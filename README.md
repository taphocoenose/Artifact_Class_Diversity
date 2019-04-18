# Artifact Class Diversity
This project contains R (www.r-project.org) and Stan (www.mc-stan.org) scripts that estimate artifact class diversity in multiple geographic regions.<br><br>

<hr>

<b>Point_Diversity.R</b><br><br>
This script requires two csv files: <i>geog_data.csv</i> and <i>point_data.csv</i>. The former csv file contains a row for every administrative unit under consideration (e.g., counties, states), with the first column consisting of unique administrative unit IDs. If the data under consideration come from the United States, I recommend using county/state FIPS codes in this column. Each succeeing column represents a geographic zone of interest, where values represent the area of the zone represented within each county. Area can be represented in any spatial unit. <i>point_data.csv</i> also begins with a column of administrative unit IDs. Every ID in <i>point_data.csv</i> should be present in <i>geog_data.csv</i> (i.e., there must be a way to tie every artifact observation to a geographic zones). Each subsequent column consists of an artifact class, with counts of that class within each administrative unit. Ensure that these two csv files and <i>point_diversity.stan</i> are in your working directory before running the script. For the examples presented by Boulanger et al. at the 2019 SAA meetings, please contact Matt Boulanger (mboulanger@smu.edu) for the appropriate csv tables.<br><br>

<i>Point_Diversity.R</i> preprocesses the csv data and proceeds to fit Stan models to <i>n</i> (user-specified) dataset permutations. For data-rich scenarios (i.e., moderate to large numbers of observed artifact counts for each geographic zone), one dataset permutation can be fitted in less than five minutes. In these scenarios, it may be possible to fit up to 100 dataset permutations on a personal computer in a reasonable amount of time. Note that this runs as a parallel process, where each model is fit on one core (the code defaults to four Hamiltonian Monte Carlo [HMC] chains per model. These are fit serially on each node). As such, quad-core machines will have significant speed advantages over dual-core machines. For larger numbers of dataets permutations (i.e., 100+), running this script on a cluster node will provide dramatic speed advantages.<br><br>

Note that model fitting speed generally depends on the amount of available data. If the dataset is dominated by regions with few observed artifacts (i.e., where most classes have 0 observations), it will take longer for Stan HMC simulations to explore the joint posterior distribution. This inefficiency will be visible in parameter estimates with low effective samples and posteriors will be vague, possibly too vague for useful inferences. In these cases, I recommend rethinking your approach, either by reducing the number of artifact classes or number of geographic zones. Although these scripts can handle uneven sampling between regions and recurring 0-count observations, it cannot provide informative inferences if there simply are not much data. For the large-scale phsyiographic example detailed by Boulanger et al. (2019 SAA meetings), we had 1435 artifact observations across 32 classes and 4 geographic zones. Most dataset permutations had at least one observed artifact for 12-32 classes in each zone. We found that this provided reasonably informative posterior estimates of class diversity in each zone.<br><br>

After the Stan model is compiled, the script will create an <i>rds</i> file in your working directory. If you would like to remove this file, wait until the R script has finished to do so.<br><br>

There are several script parameters that the user might wish to adjust depending on their research question and available computing resources. Navigate to the comment block titled 'Set main arguments' in <i>Point_Diversity.R</i>. R object <i>n_permutations</i> sets the number of dataset permutations. The default value is '6', but most users will likely wish to simulate over many more permutations than this. <i>itrtns</i> contains two integers that define the number of warmup and sampling iterations in each HMC simulation. I recommend sticking with the default 'c(5000, 5000)'. <i>adelta</i> controls the adapt_delta argument when the stan function is called. Again, I recommend leaving this value at the default '0.95'. Finally <i>div_arg</i> controls the diversity measure. The default is set to 'EH' for Shannon evenness. Other options include the Shannon Diversity Index, Simpson Diversity Index, and Inverse Simpson Diversity Index. Review the code comments in this block to choose the measure of interest.<br><br>

<i>Point_Diversity.R</i> produces two classes of <i>RData</i> files. The first of these is <i>Output_and_Diagnostics.RData</i>. This <i>RData</i> file contains the compiled Stan model (<i>Stan.model</i>), a randomly sampled example fit for one of the dataset permutations (<i>Sample_model</i>), a summary of posterior parameters for this example (<i>Stan.model.diagnostics</i>), posterior samples from this example (<i>Stan.model.samples</i>), a list of matrixes for every dataset permutation (<i>sim_point_dstrs</i>), a summary of every dataset permutation by geographic zone (<i>datasets_summary</i>), and aggregated posterior diversity samples across all model fits (<i>Post_Div_Samples</i>). Also included in this file is a dataframe  named <i>model_diagnostics</i>. This dataframe contains a summary of convergence diagnostics for every model fit, including the maximum Rhat value, minimum effective sample size (as a percentage of the total sampling iterations), maximum treedepth, number of divergences, number of parameters with Rhat values above 1.01, number of parameters with effective sample sizes below 10% of the total sampling iterations, and the number of iterations that exceeded maximum treedepth. The final column of this dataframe flags any fits with potential convergence issues with 1. For models with sparse observations (see two paragraphs up), minimum effective samples may fall below 10%, but the model has still converged. In these cases, check the other diagnostics to ensure convergence, and examine trace plots from several of the model fits.<br><br>

The second class of <i>RData</i> files produced by <i>Point_Diversity.R</i> consists of aggregated model fits (<i>Mdls_XXXX.RData</i>). These are separated from <i>Output_and_Diagnostics.RData</i> because they use a lot of memory, and in the case of many dataset permutations, may not load into an R session on a personal machine. One file is created per 100 model fits. To inspect convergence diagnostics for a specific dataset permutation fit, load the appropriate <i>Mdls_XXXX.RData</i> file into an R session.<br><br>

<hr>

<b>Plotting_Code.R</b><br><br>
This script requires that <i>Output_and_Diagnostics.RData</i> is located in your working directory. Execute <i>Point_Diversity.R</i> to generate this file.<br><br>

<i>Plotting_Code.R</i> Generates four different plots. The first plot displays posterior diversity values for up to 30 dataset permutations. Each row is a permutation, with a density polygon for each geographic region. The densities are represented by distributions fitted to posterior samples (beta distributions for measures on the (0,1) scale and gamma distributions for measures on the (0, Inf) scale).<br><br>

The second plot displays samples aggregated across all dataset permutations, with the posterior median, 50% HPDI, and 95% HPDI symbolized by white dots and bars. Due to the potentially very large number of samples across all models, the plotted samples are thinned for display (default is 1/20 samples).<br><br>

The third plot displays estimated proportions of each artifact class in each region for the example permutation model fit. In each panel, the x-axis moves across 500 posterior samples, with stacked bands representing the proportional representation of each artifact class. When posterior uncertainty is high, the bands become noisier. Band color and y-axis position are constant across the geographic zone panels in the plot. Posterior diversity estimates calculated from these proportions are plotted with points above each panel of stacked bands.<br><br>

The fourth plot displays "scree" line geometry for each geographic zone in the example model fit. Artifact classes are placed in ascending order of their posterior median proportions along the x-axis (x-axis placement changes for each artifact class across geographic zone panels). White points indicate median posterior proportions, and blue geometry illustrates uncertainty around these proportions through transparency level.
