### Variant Classifier

## Problem

In this example, we attempt to classify variants based on allele number, mode of inheritance of disease, and allele frequency.

## Dataset

The dataset is synthetic to provide a proof of concept of using machine learning to automate the assignment of weights to pieces of pathogenicity evidence.

## Origin

This example is a simplified version of one of the experiments from Bach et al.'s core PSL paper:
"Hinge-Loss Markov Random Fields and Probabilistic Soft Logic":
```
@article{bach:jmlr17,
  Author = {Bach, Stephen H. and Broecheler, Matthias and Huang, Bert and Getoor, Lise},
  Journal = {Journal of Machine Learning Research (JMLR)},
  Title = {Hinge-Loss {M}arkov Random Fields and Probabilistic Soft Logic},
  Year = {2017}
}
```

It is based on the citation-categories example from the PSL developers. https://github.com/linqs/psl-examples/tree/master/citation-categories
