### Howto merge a branch

We assume you have forked TUWien-ASC/ASC-bla into your personal repo on github, like JSchoeberl/bla. Then you have cloned this one to your computer.


First, we have to add the original git repo as a remote repo. Go to REMOTE, add it using the `+` and name it as you like (I have called it *dasOriginal*).



```{image} AddRemote.png
:width: 40%
:align: center
```

Then go to the git controls and select `branch` -> `merge branch`. Choose the branch `dasOriginal - expr`, it will merge the changes from there into your branch.

My git-graph view looks now like

```{image} GitGraph.png
:width: 40%
:align: center
```

Then, you synchonize with your upstream repo.




