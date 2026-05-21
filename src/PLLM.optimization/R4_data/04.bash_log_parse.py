import os
import pandas as pd

def parse_training_logs(log_files):
    """
    Parse training log files and extract metrics into a DataFrame.
    
    Parameters
    ----------
    log_files : list of str
        List of log file paths to parse.
    
    Returns
    -------
    pd.DataFrame
        DataFrame containing parsed metrics with columns:
        seed, model, data, batch_size, method, count, num_epochs, lr, loss, 
        epoch, type, spearman
    """
    res = pd.DataFrame()
    
    for _index, _log in enumerate(log_files, start=1):
        with open(_log, 'r') as file:
            lines = file.readlines()
            _i = 0  # Train epoch counter
            _j = 0  # Validation/Test epoch counter
            
            # Initialize variables to avoid UnboundLocalError
            _seed   = None
            _model  = None
            _data   = None
            _method = None
            _count  = None
            _epochs = None
            _lr     = None
            _loss   = None
            _batch_size = None
            
            for _line in lines:
                if " Random seed:" in _line:
                    _seed = int(_line.split(":")[-1].strip())
                if " Model:" in _line:
                    _model = str(_line.split(":")[-1].strip())
                    _model = _model.split("/")[-1]
                if " Data:" in _line:
                    _data = str(_line.split(":")[-1].strip())
                if " Method:" in _line:
                    _method = str(_line.split(":")[-1].strip())
                if " Rank count:" in _line:
                    _count = int(_line.split(":")[-1].strip())
                elif " Top layers:" in _line:
                    _count = int(_line.split(":")[-1].strip())
                if " Epochs:" in _line:
                    _epochs = int(_line.split(":")[-1].strip())
                if " lr:" in _line:
                    _lr = float(_line.split(":")[-1].strip())
                if " MSELoss:" in _line:
                    _loss = "MSE"
                elif " ListMLE:" in _line:
                    _loss = "ListMLE"
                if "Batch size" in _line:
                    _batch_size = int(_line.split(":")[-1].strip())
                    
                if "Train Spearman" in _line:
                    _train_spearman = float(_line.split(":")[-1].strip())
                    _i += 1
                    _s = pd.DataFrame(
                        data=[[_seed, _model, _data, _batch_size, _method, _count, 
                               _epochs, _lr, _loss, _i, "train", _train_spearman]], 
                        columns=["seed", "model", "data", "batch_size", "method", 
                                "count", "num_epochs", "lr", "loss", "epoch", 
                                "type", "spearman"]
                    )
                    res = pd.concat([res, _s], ignore_index=True)
                # Test data performance is corresponding to the last epoch
                if "Validation Spearman" in _line or "Test Spearman" in _line:
                    _valid_spearman = float(_line.split(":")[-1].strip())
                    _j += 1
                    _s = pd.DataFrame(
                        data=[[_seed, _model, _data, _batch_size, _method, _count, 
                               _epochs, _lr, _loss, _j, "valid", _valid_spearman]], 
                        columns=["seed", "model", "data", "batch_size", "method", 
                                "count", "num_epochs", "lr", "loss", "epoch", 
                                "type", "spearman"]
                    )
                    res = pd.concat([res, _s], ignore_index=True)
                    
        print(f"Processed {_index}/{len(log_files)} logs")
    
    return res


if __name__ == "__main__":
    # Get all log files starting with "03" and ending with "log"
    log_path = "./logs"
    R4_logs  = [os.path.join(log_path, i) for i in os.listdir(log_path) if i.startswith("03") and i.endswith(".log")]
    
    # Parse the logs
    R4_res   = parse_training_logs(R4_logs)

    # Save results (using the last lr value if available)
    if not R4_res.empty and 'lr' in R4_res.columns:
        _lr = R4_res['lr'].iloc[-1]
        R4_res.to_csv(f"04.R4.log_parse.lr{_lr}.log", index=False)
    else:
        R4_res.to_csv("04.R4.log_parse.log", index=False)
