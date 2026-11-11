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
    fig, axes = plt.subplots(3, 3, figsize=(15, 12), constrained_layout=True)
    for ax, entry in zip(axes.flat, plan['models'], strict=True):
        name = entry['name']
        source = step / 'execution/tmp/round1/fits' / name / 'run.json'
        raw = source.read_bytes()
        run = json.loads(raw)
        assert run['status'] == 'complete' and run['model'] == name
        assert run['policy'] == policy and run['policy_sha256'] == plan['policy_sha256']
        assert run['fold'] == plan['fold'] and run['seed'] == plan['seed']
        assert run['predictor_sha256'] == entry['predictor_sha256']
        sources[str(source.relative_to(step))] = hashlib.sha256(raw).hexdigest()
        for phase in ('tahoe', 'tenv'):
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
            offset = 0 if phase == 'tahoe' else policy['phases']['tahoe']['updates']
            positions = [offset + update for update in updates]
            # Plot stages separately: fine-tuning restarts from the selected Tahoe checkpoint.
            ax.plot(positions, train, 'o-', color='blue', linewidth=2)
            ax.plot(positions, valid, 'o-', color='green', linewidth=2)
            best = valid.index(min(valid))
            ax.scatter(positions[best], valid[best], s=110, facecolors='none',
                       edgecolors='black', zorder=5)
        boundary = policy['phases']['tahoe']['updates']
        ax.axvline(boundary + 150, color='#888888', linestyle='--', linewidth=1)
        ax.set_title(name, fontsize=12, weight='semibold')
        ax.set_xticks([256, 512, boundary + 300, boundary + 600],
                      ['Tahoe\n256', 'Tahoe\n512', 'tenv\n300', 'tenv\n600'])
        ax.set_ylabel('Weighted objective loss')
        ax.grid(alpha=.22)
        ax.tick_params(labelsize=9)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color='blue', marker='o', label='train'),
                        Line2D([], [], color='green', marker='o', label='valid'),
                        Line2D([], [], color='black', marker='o', markerfacecolor='none',
                               linestyle='none', label='best recorded valid (per phase)')],
               loc='upper center', bbox_to_anchor=(.5, .965), ncol=3, frameon=False)
    fig.get_layout_engine().set(rect=(0, 0, 1, .93))
    fig.suptitle('XPert Round 1 learning curves | fold 0, seed 17', fontsize=17, y=.995)
    fig.supxlabel('Update labels are phase-local; stages are not connected. Only two recorded points per phase.\n'
                  'Train: preceding-window mean; valid: checkpoint evaluation. Dashed line: phase boundary.',
                  fontsize=10)
    fig.savefig(output / 'xpert_opti_r1_learning_curves.png', dpi=180,
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
