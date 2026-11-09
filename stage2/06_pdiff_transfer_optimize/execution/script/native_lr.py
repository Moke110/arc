"""Add a fixed-step EMA snapshot to the byte-frozen Step 04 native trainer."""
from unittest.mock import patch

from artifacts import sha256, write_json
from training_inputs import json_hash

NATIVE_CHECKPOINT = '01-native-selection.pt'


def common_checkpoint(step):
    return f'02-common-step-{step}.pt'


def fit_diagnostic(native, runtime, request, recipe, directory, common_step):
    """Keep native selection immutable once it stops and continue its live fit."""
    config = native.load_recipe(runtime, recipe, request)
    if (type(common_step) is not int or not 1 <= common_step <= config['max_steps']
            or (common_step % config['validation_interval'] and common_step != config['max_steps'])):
        raise ValueError('common optimizer step must be on the frozen validation cadence')

    model = stop = None
    selection_class = native.Selection
    model_class = native.JointPerturbDiff

    class TrackModel(model_class):
        def __init__(self, *args, **kwargs):
            nonlocal model
            super().__init__(*args, **kwargs)
            model = self

    class FixedStepSelection(selection_class):
        def observe(self, step, loss):
            nonlocal stop
            if step == common_step:
                path = directory / common_checkpoint(step)
                state = {'identity': {**request['identity'], 'seed': config['seed'], 'stage': config['stage']},
                         'recipe_sha256': json_hash(config), 'best_step': step,
                         'best_validation_loss': loss, 'model_state_sha256': native.digest(model)}
                native.save(path, {**state, 'format': native.FORMAT, 'recipe': config,
                                   'model_state': model.state_dict()})
                write_json(directory / 'diagnostic.json', {**state, 'checkpoint_sha256': sha256(path)})
            if stop is not None:
                # Native still validates its EMA while driving the same fit.
                # Its selection and history end at the previously sealed stop.
                return {'step': step, 'checkpoint_refreshed': False, 'early_stop': step == common_step}
            event = super().observe(step, loss)
            if event['early_stop'] or step == config['max_steps']:
                stop = {'stop_step': step,
                        'stop_reason': 'validation_early_stop' if event['early_stop'] else 'max_steps'}
                if step < common_step:
                    return {**event, 'early_stop': False}
            return event

    def write_run(path, run):
        if path == directory / 'run.json' and stop is not None:
            stop_step = stop['stop_step']
            run['history'] = [
                {**event, 'early_stop': stop['stop_reason'] == 'validation_early_stop'}
                if event['step'] == stop_step else event
                for event in run['history'] if event['step'] <= stop_step
            ]
            seal = directory / 'native-selection.json'
            # Only seal after the native trainer has finished the stop event
            # and saved its selected checkpoint. A failure during that event
            # leaves the native selection unavailable.
            if (run['status'] == 'running' and run['history']
                    and run['history'][-1]['step'] == stop_step and not seal.exists()):
                write_json(seal, {**stop, 'monitor': run['monitor'], 'history_sha256': json_hash(run['history']),
                                  **{key: run[key] for key in ('identity', 'recipe_sha256', 'best_step',
                                      'best_validation_loss', 'model_state_sha256', 'checkpoint_sha256')}})
            if run['status'] == 'completed':
                run['stop_reason'] = stop['stop_reason']
        write_json(path, run)

    try:
        with patch.multiple(native, JointPerturbDiff=TrackModel, Selection=FixedStepSelection, write_json=write_run):
            native.fit(runtime, request, recipe, directory)
    finally:
        best = directory / 'best.pt'
        if best.exists():
            best.rename(directory / NATIVE_CHECKPOINT)
