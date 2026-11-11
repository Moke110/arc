"""Plot and audit existing Round 1 loss logs; no fitting or expression access."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Frozen Round 1 objective from model.py; avoid importing the training runtime.
LOSS_WEIGHTS = (0.2, 0.003, 0.2, 1.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step', type=Path, default=Path('stage2/09_xpert_optimize'))
    args = parser.parse_args()
    step = args.step
    plan = json.loads((step / 'execution/round1_record/plan.json').read_text())
    policy = json.loads((step / 'execution/round1_record/policy.json').read_text())
    output = step / 'output/round1'
    points, summary, sources = [], [], {}
    for entry in plan['models']:
        name = entry['name']
        source = step / 'execution/tmp/round1/fits' / name / 'run.json'
        raw = source.read_bytes()
        run = json.loads(raw)
        assert run['status'] == 'complete' and run['model'] == name
        assert run['policy'] == policy and run['policy_sha256'] == plan['policy_sha256']
        assert run['fold'] == plan['fold'] and run['seed'] == plan['seed']
        assert run['predictor_sha256'] == entry['predictor_sha256']
        sources[str(source.relative_to(step))] = hashlib.sha256(raw).hexdigest()
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
        for ax, phase in zip(axes, ('tahoe', 'tenv'), strict=True):
            result = run['phases'][phase]
            history = result['history']
            updates = [row['update'] for row in history]
            valid = [row['validation_loss'] for row in history]
            train = [sum(w * t for w, t in zip(LOSS_WEIGHTS, row['train_terms'], strict=True))
                     for row in history]
            assert len(updates) >= 2 and updates == sorted(set(updates))
            assert updates[-1] == policy['phases'][phase]['updates']
            delta = valid[-1] - valid[-2]
            at_minimum = valid[-1] <= min(valid)
            passes = delta >= 0 and not at_minimum
            summary.append(dict(model=name, phase=phase, recorded_points=len(history),
                                final_update=updates[-1], best_recorded_update=updates[valid.index(min(valid))],
                                previous_valid_loss=valid[-2], final_valid_loss=valid[-1],
                                final_change=delta, change_percent=100 * delta / valid[-2],
                                final_is_recorded_minimum=at_minimum,
                                meets_endpoint_criterion=passes))
            for row, loss in zip(history, train, strict=True):
                points.append(dict(model=name, phase=phase, update=row['update'],
                                   train_loss=loss, valid_loss=row['validation_loss']))
            ax.plot(updates, train, 'o-', color='blue', label='train', linewidth=2)
            ax.plot(updates, valid, 'o-', color='green', label='valid', linewidth=2)
            best = valid.index(min(valid))
            ax.scatter(updates[best], valid[best], s=150, facecolors='none', edgecolors='black',
                       label='best recorded valid', zorder=5)
            ax.set_title(f'{"Tahoe pretrain" if phase == "tahoe" else "tenv fine-tune"}\n'
                         f'Endpoint criterion: {"PASS" if passes else "FAIL"}; valid change {100 * delta / valid[-2]:+.2f}%')
            ax.set(xlabel='Optimizer update (within phase)', ylabel='Weighted objective loss', xticks=updates)
            ax.grid(alpha=.22)
            ax.legend(fontsize=8)
        fig.suptitle(f'{name} | Round 1 | fold 0, seed 17', fontsize=14)
        fig.supxlabel('Only two recorded points per phase. Train: preceding-window mean; valid: checkpoint evaluation.',
                      fontsize=9)
        fig.savefig(output / f'xpert_opti_r1_learning_{name}.png', dpi=160,
                    metadata={'Software': 'Stage 2 Step 09'})
        plt.close(fig)
    for filename, rows in [('learning_curve_points.csv', points), ('learning_curve_convergence.csv', summary)]:
        with (output / filename).open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    record = dict(source_sha256=sources, loss_weights=LOSS_WEIGHTS,
                  criterion='final valid >= previous valid AND final valid > minimum recorded valid',
                  limitation='Two observations per phase only; endpoint criterion is not proof of asymptotic convergence.',
                  deterministic_references='no_change and delta_mean do not train and have no learning curves')
    (step / 'execution/round1_record/learning-curves.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
