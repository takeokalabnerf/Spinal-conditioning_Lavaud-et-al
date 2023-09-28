# Two distinct inhibitory neuronal classes govern acquisition and recall of spinal sensorimotor learning - Figure replication

## Extraction of the data
Data are extracted from the curated dataset (provided on demand, for further details, see "Preprocessing of neural recordings" in the Methods section). Selected information about single-unit are stored in a matlab table (frM or frY), for the Learner or Control animals respectively.
The following code allow the creation of the tables:
```sh
Creation_table_Learners_Controls.m
```
Each table contains the following information:
- Recording: Recording number
- Id_horr: Cluster number of the unit given by kilosort.
For the following, please note that the “early phase” and “middle phase” are cuts of the “pre-criterion phase” that are not used in the paper. The “late phase” correspond to the “post criterion phase” (see the methods section: "Single-unit analysis of genetically undefined spinal neurons").
- Tbin: Time bin representing for each row, the start and the end of the baseline (or also called in the codes spontaneous phase, the pre-criterion (early phase + the middle phase) and the post-criterion (late phase). 
- Fr_spont: 1-second bin of the firing rate during the 10 minutes of the baseline/spontaneous period.
- Fr_average_spont: Average of Fr_spont.
- Fr_early: 1-second bin of the firing rate during the early phase.
- Fr_middle: 1-second bin of the firing rate during the middle phase.
- Fr_late: 1-second bin of the firing rate during the late phase.
- Fr_after: 1-second bin of the firing rate during the 10 minutes following the end of the conditioning.
- Fr_average_after: Average of Fr_after.
For the following, z-score are calculated as
> z = (x-μ)/σ
where x is the raw value of the firing rate, μ is the average of the firing rate during the spontaneous period, and σ is the standard deviation of the firing rate during the spontaneous period
- zscore_early_full: Zscore of the early phase.
- zscore_mid_full: Zscore of the middle phase.
- zscore_late_full: Zscore of the late phase.
  For the following, once z-score are calculated, they are getting assigned a value : 1 or 10 if zscore>2 for 50 or 20% or the phase respectively, 2 or 20 if zscore<-2 and 0 zscore<2 and zscore>-2. . For detailled calculation, see the methods section: "Single-unit analysis of genetically undefined spinal neurons"
- zscore_group_EandMmerge: Zscore assignement score for the early and middle phase merged.

Not used in our study :
- zscore_EandMmerge_up: value of the average zscore increase.
- zscore_EandMmerge_down: value of the average zscore decrease.
- zscore_group_lateZscore assignement score for the late phase.
- zscore_late_up: value of the average zscore increase.
- zscore_late_down: value of the average zscore decrease.

We used the assigned value based on the zscore to categorize these units as either "upregulated", "downregulated," or "no change" based on this z-score value, with a "Category" number. 1 and 2 are an increase during both the early+middle and the late phase, 3 and 4 are an increase during the early+middle phase only, 5 and 6 are an increase during the late phase only, 7 is a decrease during the early+middle phase and the late phase, 8 is a decrease during the early+middle phase only, 9 is a decrease during the late phase, 10 are units that are not changing their activity during any phase, and finally 0 are the units we excluded (mix on increase and decrease <2% of dataset).
- Category: category assignement of the single unit (see above).
- Depth: estimated depth of the single units.
- Optotagged: binary number, 1 if the unit is triggered during photo-activation of Ptf1aON neurons (see "
Opto-tagging and analysis" in the Methods for more details).
- Strong_firing: Neurons with low firing rate across the trial < 0.05 Hz (binary number = 0) were considered non changing. 

Table frM and frY are extraction of single units coming from Ptf1aON animal (recording number >10 for learner and >9 for controls). frM_En1 and frY_En1 defined En1ON animals. frM_vGat and frY_vGat defined VgatON animal. frM_vGlut and frY_vGlut defined vGlutON animals. Finally, frM_Ptf1a_CNO defined Ptf1aON mice that have been injected with CNO during recording (see "Pharmacogenetic manipulations" in the Methods for more details).

For second order units (see "Categorization of input afferent types" in the Methods for more details), other tables are created call frM_UOI, frY_UOI (for Learner and Controls Ptf1aON animals respectively) and frM_UOI_CNO (for Learner Ptf1aON animal recorded with injection of CNO, see "Pharmacogenetic manipulations" in the Methods for more details).
The following code allow the creation of the tables:
```sh
Second_Order_Units.m
```
_UOI table are composed of different parameters:
- Recording: Recording number
- Id_horr: Identification number of the single unit 
- Depth: estimated depth of the single units.
- Category: category assignement of the single unit.
For the following, please see "Categorization of input afferent types" and "Determination of response probability to conditioning cues" in the methods for more details.
- latency_all: returned the latency of each responses peak
- Ptf1aON: returned 1 if the units is a Ptf1aON unit.
- vGlutON: returned 1 if the units is a vGlutON unit.
- vGatON: returned 1 if the units is a vGatON unit.
- SetOfSpikeStart: % of spikes per stimulation for the first 10% of the cues
- SetOfSpikeEnd: % of spikes per stimulation for the last 10% of the cues
- ReliabilityDelta: difference between SetOfSpikeEnd and SetOfSpikeStart

## Figure 3

### 3C
For the raw spike detection see Sans-Dublanc, Chrzanowska et al., 2021 (Sans-Dublanc A*, Chrzanowska A*, Reinhard K, Lemmon D, Nuttin B, Lambert T, Montaldo G, Urban A, Farrow, K. Optogenetic fUSI for brain-wide mapping of neural activity mediating collicular-dependent behaviors. Neuron. 2021 Apr 22;S0896-6273(21)00238-5. doi: 10.1016/j.neuron.2021.04.008.).

### 3D 
For the spike per second and events figure, the code is coming from Cortex Lab (https://github.com/cortex-lab/phy).

### 3E
For this figure, please run:
```sh
Category_of_activity_Ptf1a.m
```
It will extract from the right table, the depth of only the Ptf1aON units.

### 3H
For this figure, please run:
```sh
Analysis_Zscore_Ptf1a.m
```
It will plot the 10 minutes zscore of Learner and Controls and display the difference in firing rate.

## Figure 4

### 4B, 4E
For this figure, please run:
```sh
Second_Order_Units.m
```
It will allow you to have different example of units with different reponse latencies.

### 4D, 4F
For this figure, please run:
```sh
Analysis_Reliability_Ptf1a_CNO.m
```
It will return the reliability of only Ad/C of the Ptf1aON mice recorded with injection of CNO (see "Pharmacogenetic manipulations" in the Methods for more details).

### 4F, 4H
For this figure, please run:
```sh
AdC_P_responsive_units.m
```
It will focus on the analysis of only Ad/C responses or Ad/C contigent with P responses. It will also returned the extent in reliability change.

## Figure supplementary 5

## Figure supplementary 6

### S6B

### S6C, S6D
For this figure, please run:
```sh
Change_Zscore.m
```
It will give you the proportion of units with upregulated, downregulated or no change in their activity, and their depths.

## Figure supplementary 7

### S7A
For this figure, please run:
```sh
Second_Order_Units.m
```
It will allow you to have different example of units with different reponse latencies.

### S7B, S7D
For this figure, please run:
```sh
depth_latency_tables.m
```
It will give you the depths of differents units based on their response latency.

### S7C
For this figure, please run:
```sh
second_order_type_multi.m
```
It will give you the proportion of units that have only one response or multiple responses.

![image](https://github.com/takeokalabnerf/Spinal-conditioning_Lavaud-et-al/assets/143432651/dbfefdbe-28be-46d3-9fdf-a399b6766b95)
