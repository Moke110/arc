"""Audit the cold-drug-dose estimand against its frozen split ledger."""
import math

IDENTITIES = ['profile_id', 'drug_parent', 'dose_um']
PARTITIONS = ('train', 'validation', 'test', 'excluded-held-parent', 'excluded-held-dose')


def audit_coverage(ledger, scored_by_fold):
    """Verify exact test coverage and return three fold-level audit records."""
    auxiliary = ledger.loc[ledger.scenario.eq('cold-drug-dose')]
    expected_folds = {0, 1, 2}
    if set(auxiliary.fold) != expected_folds or set(scored_by_fold) != expected_folds:
        raise ValueError('cold-drug-dose requires exactly three outer folds')

    population = auxiliary.loc[auxiliary.fold.eq(0), IDENTITIES].reset_index(drop=True)
    if population.profile_id.duplicated().any():
        raise ValueError('cold-drug-dose population has duplicate profile identities')
    doses = sorted(population.dose_um.unique())
    if len(doses) != 3 or any(not math.isfinite(dose) or dose <= 0 for dose in doses):
        raise ValueError('cold-drug-dose requires three positive finite source doses')

    held_groups = set()
    records = []
    for fold, held_dose in enumerate(doses):
        rows = auxiliary.loc[auxiliary.fold.eq(fold)].reset_index(drop=True)
        if not rows[IDENTITIES].equals(population):
            raise ValueError(f'cold-drug-dose fold {fold} has changed frozen population')
        if set(rows.partition) != set(PARTITIONS):
            raise ValueError(f'cold-drug-dose fold {fold} is missing a partition')
        expected_test = rows.partition.eq('test')
        held_parents = set(rows.loc[expected_test | rows.partition.eq('excluded-held-parent'), 'drug_parent'])
        if held_groups & held_parents:
            raise ValueError('cold-drug-dose held parent groups overlap')
        held_groups.update(held_parents)
        parent_held = rows.drug_parent.isin(held_parents)
        dose_held = rows.dose_um.eq(held_dose)
        masks = {
            'test': parent_held & dose_held,
            'excluded-held-parent': parent_held & ~dose_held,
            'excluded-held-dose': ~parent_held & dose_held,
        }
        if any(not rows.partition.eq(partition).equals(mask) for partition, mask in masks.items()):
            raise ValueError(f'cold-drug-dose fold {fold} has incorrect test or excluded strips')
        train_parents = set(rows.loc[rows.partition.eq('train'), 'drug_parent'])
        valid_parents = set(rows.loc[rows.partition.eq('validation'), 'drug_parent'])
        if (train_parents & valid_parents
                or len(valid_parents) != math.ceil(0.1 * len(train_parents | valid_parents))):
            raise ValueError(f'cold-drug-dose fold {fold} has incorrect validation parent coverage')

        scored = scored_by_fold[fold][IDENTITIES].reset_index(drop=True)
        expected = rows.loc[expected_test, IDENTITIES].reset_index(drop=True)
        if not scored.equals(expected):
            raise ValueError(f'cold-drug-dose fold {fold} scored identities differ from frozen test rows')
        n_test = int(expected_test.sum())
        records.append(dict(
            fold=fold, held_dose_um=float(held_dose), n_held_parents=len(held_parents),
            n_population_samples=len(population), n_test_samples=n_test,
            n_test_parents=int(rows.loc[expected_test, 'drug_parent'].nunique()),
            n_train_samples=int(rows.partition.eq('train').sum()),
            n_validation_samples=int(rows.partition.eq('validation').sum()),
            n_excluded_held_parent=int(masks['excluded-held-parent'].sum()),
            n_excluded_held_dose=int(masks['excluded-held-dose'].sum()),
            test_population_fraction=n_test / len(population),
            test_scope='held-parent × held-dose intersections only',
        ))
    if held_groups != set(population.drug_parent):
        raise ValueError('cold-drug-dose held groups do not cover the source parents')
    return records
