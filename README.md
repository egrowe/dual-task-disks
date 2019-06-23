# dual-task-disks
Dual-task similarity ratings experiment using red/green disks

We used a dual-task paradigm to compare performance and perceptual similarity ratings during the central and peripheral single-tasks independently and when both were performed simultaneously. In all blocks, the visual stimuli remained the same, only the task instructions and task-relevance of the stimuli changed.

The peripheral stimulus set consisted of eight red/green bisected disks that varied in their degree of rotation from 0o to 315o in 45o increments. Overall, there were 64 possible disk pairs. In some blocks of the Main Experiment only a subset (2) of these disks was presented to participants (i.e. only the vertically bisected red/green and vertically bisected green/red disks). In each trial, two of these 8 disks were randomly selected and presented in one of two diametrically opposing positions in the periphery of the screen. The disks took up 2.5o visual angle.The two possible locations that the disk pairs could appear was at the corners of an imaginary 8o x 10o rectangle around the centre of the screen. After a short temporal delay, each of the red/green disks was replaced with an equally sized disk containing a randomly generated Mondrian pattern. The time delay between the stimulus pairs and the mask was determined via the QUEST staircase procedure.

We used a QUEST staircase procedure to determine the stimulus onset asynchronies (SOAs) for central and peripheral targets and masks. Using the QUEST (Quick Estimate of Threshold, Watson & Pelli, 1983) staircase procedure, we were able to reduce this training time by robustly adjusting the SOAs in a trial-by-trial manner for the central and peripheral (SOA setup) single-tasks according to each participant’s abilities (set at 70% accuracy in both of the single-tasks).

See the attached figures for an in-depth explanaton of the block structure, stimulus sets and scripts used to run these parts of the experiment.


NOTE: This script also contains commented out code for including different 'offsets' of the disks, from half-half and full-offet. To include these conditions, uncomment these parts in ExpParameters scripts. The variable possStim then becomes 1:24 (all 24 possible disk stimuli)
