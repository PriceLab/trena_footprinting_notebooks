# trena_footprinting_notebooks
notebooks and code for ensemble learning of footprint classification

## instructions

There are a handul of central notebooks for manipulating data and running models, and some auxillary notebooks for explatory analysis.

### Main notebooks (and the order in which to run them)

- construct_motif_dataset.ipynb: sample fimo motifs, a certain percentage of which overlap chipseq hits
- merge_motif_dataset_hint_wellington.ipynb: merge in hint and wellington scores to the above data set
- annotate_merged_dataset.ipynb.ipynb: add things like gc content, tss distance, TF class to above data set
- whole_genome_classifier_footprint_hit_subset: train decision tree classifier on the above data set, but only on instances where there is a nonzero hint, wellington, or chipseq score
- run_saved_model.ipynb: example of how to load a pre-trained model and use it for classification

### auxillary notebooks

- evaluate_motif_dataset.ipynb: check out some stats on the inital data set construction
- per_TF_classifier.ipynb: train models on TFs individually. aggregate performace is worse than training all at once.
- whole_genome_classifier_naive_set.ipynb: same as whole_genome_classifier_footprint_hit_subset.ipynb above, except don't filter out instances where there is no hint, wellington, or chipseq hit.  trees show less advantage here, though the problem is more artificial.
- old_notebooks: nothing of interest, just keeping around until i'm sure i can delete with impunity

