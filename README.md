# REINVENT benchmarking

Code is adapted from [https://github.com/MarcusOlivecrona/REINVENT](MarcusOlivecrona/REINVENT).

Updated for use with PyTorch v1.10.0. Made for use on CC cluster.


## Environment

All required imports and modules are found in `make_env.sh`. Create virtual environment using 
```
./make_env.sh
source ~/env/reinvent/bin/activate
```

## Prepare data

Put data into the `data/` directory. Convert the smiles of data into `data/data.smi` by running
```
python prepare_data.py
```

You may need to modify the script depending on the format of dataset.

## Train model

Pretrain the model for specified `num-epochs` by running:
```
python train_prior.py --num-epochs 50
```

You can view progress bar by toggling flag `--verbose`. No early stopping is implemented. Learning is decreased every 50 steps.

## Fitness function

The fitness function is a class in `custom.py`. Modify the fitness function by changing the `fitness_function` function, which takes in a SMILES and outputs a float value. Return `None` for invalid smiles.

The values will be automatically scaled by a sigmoid. For sampled molecules *m* with fitness *f*, the score is (invalid molecules have 0 reward and do not affect the agent)

![](extras/scoring_function.png)

where *a* is the average of all known fitnesses from previous generations, and *b* is the slope of the sigmoid, calculated from a threshold. 

For a set *F* of fitnesses from previous generations

![](extras/parameters.png)

where *c* is the specified threshold (default 0.8). This means for the current fitnesses *F*, the maximum fitness will map to 0.8 on the sigmoid. The smaller the threshold *c*, the stronger the reward for fitnesses larger than *F*.

All other parameters are left as default from REINVENT 1.0.

## Run REINVENT

Run the REINVENT for generations = `num-steps` with evaluations = `batch_size` per generation. Run the script
```
python main.py --num-steps 10 --batch-size 500
```
