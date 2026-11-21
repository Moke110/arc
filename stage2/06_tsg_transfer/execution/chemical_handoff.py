"""Transfer only the native model's fingerprint embedding module."""


def initialize_chemical(source, recipient):
    """Strictly load all embedding parameters and buffers, including BN state."""
    state = source.feat_embeddings.state_dict()
    target = recipient.feat_embeddings.state_dict()
    if not state or set(state) != set(target) or any(
            state[name].dtype != target[name].dtype or state[name].shape != target[name].shape
            for name in state):
        raise ValueError('incompatible chemical modules')
    recipient.feat_embeddings.load_state_dict(state, strict=True)
    return tuple('feat_embeddings.' + name for name in state)
