def __getattr__(name):
    if name == 'XPertNet':
        from src.model_xpert import XPertNet
        return XPertNet
    elif name == 'MyDataset':
        from src.dataset import MyDataset
        return MyDataset
    elif name in ('get_metrics_new', 'get_metrics_infer'):
        from src.metrics import get_metrics_new, get_metrics_infer
        return get_metrics_new if name == 'get_metrics_new' else get_metrics_infer
    elif name == 'get_metrics_distribution':
        from src.eval_metrics import get_metrics_distribution
        return get_metrics_distribution
    elif name in ('EarlyStopping', 'set_random_seed'):
        from src.utils import EarlyStopping, set_random_seed
        return EarlyStopping if name == 'EarlyStopping' else set_random_seed
    raise AttributeError(f"module 'src' has no attribute {name}")
