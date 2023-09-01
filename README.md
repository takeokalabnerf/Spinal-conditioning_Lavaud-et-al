<h1> Two distinct inhibitory neuronal classes govern acquisition and recall of spinal sensorimotor learning </h1>
<h2> Code repository for kinematic analysis </h2>
<h3> Step 1: Data </h3>
All the data you need to run this cord are stored under:
> ConditioningExample/DataExp1/Day1/pair1
It is made of three raw kinematic files returning the position of 5 joints (Foot, ankle, knee, hip and crest) along the three axis X, Y and Z
<h3> Step 2: Format the kinematic component of each joint and for Learner and Control </h3>
In order to separate Learner and Control kinematics and to remove the potential artefact from the dataset, run:
> sorting_raw_kinematic_file.m
It will save the formatted kinematics in "ConditioningExample/DataExp1/AnalyseDay1/pair1".
