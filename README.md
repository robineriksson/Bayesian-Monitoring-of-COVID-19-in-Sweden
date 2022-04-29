# Bayesian-Monitoring-of-COVID-19-in-Sweden

## About
This README describes the code we use in the paper with the same name
and allows for a replication of the results in the paper.

## Authors
* Robin Marin (robin.marin 'at' it.uu.se), Primary contact
* Håkan Runvik (hakan.runvik 'at' it.uu.se),
* Alexander Medvedev (alexander 'at' it.uu.se), and
* Stefan Engblom (stefane 'at' it.uu.se), Secondary contact

# Installation
1. Clone the repository
2. Make sure that the dependencies are installed and initialized
3. Initialize by starting Matlab and running `startup` in the main
   directory

# Posteriors on file
The posterior files that are supplied in `/inference/results/KLAM/`
have all been reduced in size to 100 parameter samples for storage
size reasons. See Section *KLAM* on how to generate a full size
posterior and reproduce the results in the paper to an equivalent
precision.

# Code
Basic functionallity to generate the results in the paper. If further
explanation is wanted, see the extensive Matlab `help` for each
function. See the *Examples* section below for explicit usecases and
figure generation.

## Data (/data/)
To access data, run `d = loadData(rep)` where `rep` is a specificed
source, see `help loadData`. *Note:* not all sources listed are
distributed. Raw `.csv` files are stored under `/data/sources/REP` for
the source `REP`. We often use pre-processing steps on the data before
performing any calculations, e.g.,
```
Data = loadData('C19');
Data = polishData(Data,'D','Dinc',1);
Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});
```

## Kalman filter (/kalman/)
The main function for the Kalman filter is `C19filt(...)`. The
function makes further calls to subordinate functions found in
`/kalman/auxil/` that handle the construction and the recursive
propagation of the filter.

An evaluation of the prediction capabilities of the posterior Kalman
filter is obtained by `C19filt_lag(...)`, which is the Kalman filter
extended to enable `k`-day ahead ("lag") predictions, i.e., without data
observations.

## Prior (/inference/)
`rates = priorenger(N)` returns `N` samples from the prior, and `dens
= priorenger(rates)` computes the prior density of the samples in
`rates`.  The marginal priors are constructed in the `genprior`
script.

## Posterior (/inference/)
The posterior is stored and accessed using `rates =
posteriorenger(file)` where `file` is a posterior file.  We have
stored our generated posterior files under `inference/results/`. The
function includes a functionallity of referencing a list of posterior
files and sampling from a weighted average of them, e.g., the Swedish
national average weighted by regional population size.

## KLAM (/inference/KLAM/)
We generate samples from the approximate posterior by an
implementation of Kalman Likelihood Adaptive Metropolis (KLAM). We
supply the wrapper script `klam_parachain` which calls the underlying
functions (`klam_init`, `AM`, `logKL`) using Matlab's `parfor`
framework. Generating samples are expensive and we include the code
for reference, but the results we have generated can also be accessed by
the posterior functionality discussed above.
```
nMCMC = 1e3;        % number of samples per chain
savetofile = false; % save the posterior or not
                    % (in either case RATES is created)
reg = 2;            % Uppsala region
gelmanRub = true;   % convergence score
burnin = 1e0;       % burnin length to remove from chain
runs = 4;           % number of parallel chains
register = 'C19';   % source of data
klam_parachain;

% reproduce the results in the paper
nMCMC = 5e4;
burnin = 1e3;
reg = 1:21;
savetofile = true;
klam_parachain;
```

## DynOpt (/inference/dynamic_beta/)
The temporal upscaling ("bootstrap on the margin") of the 4-week
reproduction number is generated by `dynamic_beta_ML`. Note that if
`savetofile = true`, the script might overwrite files.
```
reg = 2;          % Uppsala
evalplot = true;  % diagnostic plots
savetofile = false;
register = 'C19'; % source of data
dynamic_beta_ML;  % generates R_t for Uppsala upscale,
                  % see the variable R_POST

reg = 1:21 % all regions
evalplot = false;
% generates the files for all the regions:
dynamic_beta_ML;

% generates the weighted national average posterior:
dynamic_beta_ML_all;
```

## Weekly predictions (/weekly/)
We perform k-day ahead prediction using the posterior Kalman filter in
`weekly_prediction`, for an example, see Fig. 3 in the paper.

## Bootstrap (/URDME/)
The boostrap replicates we use for assessing robustness are generated
using URDME. The URDME-model script `covid19enger` defines the rates
and the model itself and `covid19enger_run_post` simulates the
samples. Our posterior samples are generated using a national weighted
mean posterior. See the example for Fig. S8 below for some sample
code. The following example will generate the bootstrap replicate
samples, infer one posterior using KLAM on that data for the Uppsala
region, and finally perform the upscaling of the reproduction number
on that posterior. *Note:* these results take a long time to generate.
```
savetofile = true;
dynamic_beta_ML_all; % if not already run
regplot = [];
URDMEsampling; % generates data

reg = 2;
register = 'URDME1'
klam_parachain
urdme = 1;
dynamic_beta_ML;
```

# Examples (paper/predict/)
All examples below are associated with a display item in the paper.

Prior to generating the figure, simulation files are needed for some
of the figures/examples listed below. In order to run them properly we
suggest you first generate them by first running the `laggen` script.
It will take a while, but it needs only to be called once. The
`laggen` call should be called as follows
```
% suggested approach to reproduce the results
type = 1;
reg = [1 2 22];
laggen;
type = 2;
reg = [2]
laggen;

% if all data is wanted
type = 1;
reg = 1:22;
laggen; % filter output + 14 days of prediciton

type = 2;
laggen % Continous 7-day ahead prediction
```

## Figures (paper/predict/img/scripts)
### Fig. 2:
Illustrate the marginal prior and posterior with
`prior_posterior`. The illustration can be the weighted national
average or any combination of regional posteriors, see the variable
`reg` which include all if set to `[1:21]`, or only Uppsala if `[2]`.
```
savetofile = false;
reg = [1:21];
prior_posterior; % generates the weighted national average
                 % (as in the paper)

reg = [2]
prior_posterior; % Uppsala only
```

### Fig. 3:
7-day ahead predictions per region, and exemplified for Uppsala in
`lagplot`. The function unpacks prediction samples generated using
the `weekly_predict` wrapper `laggen`

```
savetofile = false;
generateData = 0; % if the laggen call was already made
reg = [2];        % Uppsala figure, [1] for Stockholm and so on
lagplot;
```

### Fig. 4:
The proportion of recovered invididuals per region, and Stockholm in
particular, see `recovered`. Stockholm in particular because we have
included validating sources which only consider that region. The
validating sources are included as .csv: [`fhmAnti`, `fhmAntiGivare`,
`RecPaperEstimate`], see the paper for the sources.  The call is
simply
```
savetofile = false;
recovered; % generates the plot in the paper
```

### Fig. 5 (& 7)
The posterior reproduction number (4-week constant) and the marginal
boostrap (daily) per region can be illustrated by `Rposterior`. The
region of interest is specified by `reglist`. In Fig. 5 `= [2]`, and
in Fig. S2 `= [1 10 12 8 9 19]`. The 4-week posterior is illustrated
using a boxplot, the daily as a red line, and for comparison we also
include the estimate given by PHA per region.
```
clear reg
savetofile = false;
Rposterior; % generates the figures for the 7 selected
            % regions in the paper

reg = [2];
Rposterior; % only generates the figure for Uppsala
```

### Fig. 6:
As the posterior Kalman filter can give an estimate of the proportion
of recovered individuals, similarly the filter can also give an
estimate of the number of symptomatic and/or the symptomatic
incidence. The latter is often what screening studies via testing try
to uncover. In `IincFac` we do exactly this, we query the model for
the symptomatic incidence and compare it to the positive tests by PHA
in Stockholm and Uppsala. Figures can be generates in batch, or by
region.
```
savetofile = false;
reg = [1] % Stockolm
IincFac;

% (the figure names are re-used and we suggest saving:)
savetofile = true;
reg = [1 2] % Stockholm and Uppsala
IincFac;
```

### Fig. 9:
The pre-processing of the data (as discussed under Results → Data) and
further detailed in the SI smoothes the distribution of deceased
incidence per weekday to achieve one that is closer to uniformly
distributed. The illustration in the SI is reproducable in
`weekday_smoothing`.
```
savetofile = false;
reg = 1; % Stockholm
weekday_smoothing;

reg = 2; % Uppsala
weekday_smoothing;
```

### Fig. 11:
Our exploration of the prior predictive distribution is given by
`priorpred`. The parameters are here direct samples from the prior in
the same way as the posterior is used in `laggen`. A first inital run
should include `gendata = true` and then it can be set to `false`. The
first run generates the prior sample predictions and save the file (it
can be large, therefore we do not include it). After it has been
generated, the file can simply be loaded and the predictive
distribution can be explored further.
```
savetofile = false;
regen = 1;    % generate the samples
reg = 2;      % Uppsala
Nprior = 1e3; % number of prior samples
priorpred;

regen = 0;
priorpred;    % same figure format, but loading the prior samples from file
```

### Fig. 12:
The daily estimate of β is expensive to run; superlinear complexity in
the number of days *K*. We illustrate the splitting of horizions to
make the calculations feasable in `HorizonSplitCompare`. *Warning:*
this script takes a long time.
```
savetofile = false;
reg = 2; % Uppsala region, as in paper
HorizonSplitCompare;
```

### Fig. 13:
The posterior for all regions is tricker to visualize all at once. In
`errorbarsRegion` we give the posterior mean ± 1 std per region. This
gives a quick overview of potential outlier regions but also displays
the nice agreement between regional results.
```
% load data into memory and generate figure
clear ratenames
savetofile = false;
errorbarsRegion;

% alternatively, generate figure from memory:
errorbarsRegion;
```

### Fig. 14:
The bootstrap samples are generated using URDME. The function
`URDMEsampling` generates the samples for all regions, stores the
file, and generates the figures. The latter allows for visual
comparison with the underlying data. There will be some mismatch since
the posterior we use in the simulations is a national average (except
for the reproduction number) and therefore, if the regional posterior
differed significantly from the weighted national average, we then
expect the simulations to differ somewhat. *Note:* URDME needs to be
modified a bit to work on the Windows platform, meaning that this
example is not possible to quickly generate on that system.
```
savetofile = false;
reg = [1:21];  % national average
regplot = [2];
URDMEsampling; % only the Uppsala figure

% all the figures in the paper:
clear regplot
URDMEsampling;
```

### Fig. 15:
The reproduced posterior samples from the bootstrap replicate data can
be visualized together with the weighted national average (cf. Fig. 2)
as another marginal density in `prior_posteriorURDME`.
```
savetofile = false;
reg = [1:21];         % national average
prior_posteriorURDME; % same figure format as in the paper

reg = [2];
prior_posteriorURDME; % only for Uppsala.
```

### Fig. 16:
Similarly to how the boostrap replicate posterior was compared to the
marginal posterior and prior in Fig. S9, we can include the replicate
mean for the reproduction number in `RposteriorURDME`. We only include
the files necessary for the Uppsala region plot (as given in the
paper).
```
savetofile = false;
RposteriorURDME;
```

### Fig. 17:
To compare the posterior Kalman filter predictor with something more
basic, we construct an Autoregressive (AR) model in `arx_fit`. The AR
model considers data `(H,W,D)` in an expanding window and makes 7-day
ahead predictions. The performance is then evaluated on the frequency
of "inside [X%] CrI" and the NRMSE per the days we have recorded for
the posterior filter. The reasoning for only doing it on the recorded
dates is that this ensures that the Kalman filter had not seen any
future data when making predictions. The code also generates Tab. S5.
```
savetofile = false;
reg = 2; % Uppsala region
arx_fit; % first time loads solution into memory

arx_fit; % no new fitting, only plotting
clear ypred;
reg = 1; % Stockholm data
arx_fit; % fit anew.
```

## Tables (paper/predict/tab/scripts)
### Tab. 1
We describe in the paper how the techniques from the paper was used
for weekly prediction in reports communicated with the local
authorities. The results from those reports are summarized in a table
and `weekly_eval` extracts the results.
```
savetofile = false;
weekly_eval; % reads the table
tableWeekly  % the final LaTeX-table
```

### Tab. 2
Like the reproduction number, the Infection Fatality Rate (IFR) is a
4-week dynamical parameter. We present the IFR using the
`IFRtableURDME` script. The IFR estimator per 4 week is, however,
somewhat noisy and we decided to report bi-monthly estimate instead. A
downward trend could be seen quite early on. We give our estimate
(posterior sample) for Stockholm and Uppsala and also the weighted
national average. The posterior CrI is computed as marginal
quantiles. The table include footnotemarks (§ and ‡) and we use these
marks to indicate that the estimates should be considered less robust
due to a somewhat worse bootstrap match. See the discussion in
Materials & Methods. *Note:* we only distribute a thinned
posterior due to data limitation. To fully reproduce the results,
larger sample size is needed.
```
% generates a table of the same format as in the paper
clear reg bimonthly
verb = 1;
reg = [1 2 22]
IFRtableURDME; % Stockholm, Uppsala, Sweden

reg = [1:22];  % all regions and Sweden
IFRtableURDME

reg = [2];
bimonthly = false;
IFRtableURDME; % Uppsala, but IFR per month
```

### Tab. S4
With the estimated bias we compute some uncertainty statistics:
Coefficient of Variation (CoV), Coefficient of Bias (CoB), and
normalized root mean square error (NRMSE), per region in
`bootstraptable`. We estimate the bias for the regional posteriors by
using the known bias of the bootstrap replicates. The bias is computed
per parameter, and is intended as a basic measure of the bias for the
entire region (over all parameters). Some deeper investigation into
the working of the table construction can be achieved by examining
`posterror` which processes the posterior files that `bootstraptable`
supplies it with.
```
saveetofile = false;
bootstraptable; % generates the table in the paper.
```

## Other
### CFR conditioned on [I,H,W]
In the paper we give estimates for the CFR conditioned on I, H, and
W. These esimates are generated using `spectral`.
```
reg = 1:21; % national average
sig = 3;    % # of significant digits in rounded result
spectral;

reg = 2;    % Uppsala CFR
spectral;
```

### Death per source compartment
We include also the script for estimates of the number of death per
source compartment `(I,H,W)`: `Deadsplit`. In order to run the script
a simulation output generated by an extended version of our posterior
model that keeps track of the different sources that the diseased came
from. If necessary, we will share the the additional code
needed. Please contact the owner of the repository to make that
inquiry. The current repository includes one such simulation file, but
the code to generate the simulation file is not.
```
Deadsplit; % run spectral and unpacks the simulation file
```

# Dependencies
* stenglib: https://github.com/stefanengblom/stenglib
* URDME: https://github.com/URDME/urdme
  (Windows platform not currently supported)
* MATLAB (>= release 2019a)

## Tested on
* Linux (Pop!_OS 20.10, 64 bit), MATLAB 2021a
* Linux (Pop!_OS 21.10, 64 bit), MATLAB 2021a
* macOS Monterey Version 12.3, MATLAB 2021b
* Windows 10 Education 64 bit, MATLAB R2019a
