"""All-dose no-dose adapter: unchanged native TranSiGen arithmetic and RNG."""
from . import shared
from tenv5k.model_transigen import Adapter, CONFIG as NATIVE_CONFIG, SOURCE_DIGEST, configure_backend


def adapter(source,method='no-dose',transform=None):
    if method!='no-dose' or transform is not None:raise ValueError('unsupported method or transform')
    return Adapter(source)


def fit_transform(method,doses):
    if method!='no-dose':raise ValueError('unsupported dose method')
    return None


def configuration(method,transform):
    if method!='no-dose' or transform is not None:raise ValueError('unsupported dose method')
    return dict(NATIVE_CONFIG,method=method,method_identity='tsg-dose-no-dose-native-f32-v1',transform=transform,
                initialization='native selective Xavier; unchanged no-dose constructor order')
