import types
import torch
import torch.nn as nn


def device_info():
    print(f"CPU is available: {torch.cpu.is_available()}")
    print(f"CUDA is available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        for i in range(torch.cuda.device_count()):
            print(f"CUDA device at index {i} is named {torch.cuda.get_device_name(i)}")


def log_likelihood_loss(input, *args):
    """
    Return the negative log-likelihood loss of the input, which is a tensor of
    log-likelihoods. This method does not use target values.
    """
    return -input.sum()


def which_loss(loss):
    """
    Return the specified loss function (one of "mse" , "kl_div", or "logl"). If the
    input loss is a function type, then it is returned as is.
    """
    fn_types = (
        types.BuiltinFunctionType,
        types.BuiltinMethodType,
        types.FunctionType,
        types.LambdaType,
        types.MethodType,
    )
    if hasattr(loss, "forward") or type(loss) in fn_types:
        return loss
    if type(loss) == str:
        return {
            "mse": nn.MSELoss(),
            "kl_div": nn.KLDivLoss(reduction="batchmean", log_target=True),
            "logl": log_likelihood_loss,
        }[loss]
    raise NotImplementedError(f"nothing yet for {loss}")


def random_dist(n_matrices):
    """
    Return a random vector of length n_matrix with entries in [0,1] that sum to 1.
    """
    rand = torch.rand(n_matrices)
    return rand / rand.sum()
