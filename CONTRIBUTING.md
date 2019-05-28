### Contributing to stripy

We welcome contributions to stripy, large or small [*](#footnote1). That can be in the form of new code, improvements to the documentation, 
helping with a missing test, or it may even just be pointing out a bug or potential improvement to the team. 

For bugs and suggestions, the most effective way to reach the team is by raising an issue on the github issue tracker. Github
allows you to classify your issues so that we know if it is a bug report, feature request or feedback to the authors.

If you wish to contribute some changes to the code then you should submit a *pull request*. On github we can review the code that you are contributing and discuss it before the changes are merged into the development version of the code. Before embarking on changes to the code, please first take a look at the [`dev` branch](https://github.com/underworldcode/strip/tree/dev) which is where unreleased changes are staged: we may have 
been working on something similar already. 

#### How to create a Pull Request (PR)

Create your own fork of the project on github: log into your github account, navigate to the [stripy repository](https://github.com/underworldcode/strip/tree/dev)
and click on the *fork* button at the top right corner of this repository. You can find detailed instructions and a discussion of how to fork a repository, clone it locally and work on the changes
[in the GitHub guides](https://guides.github.com/activities/forking/). 

#### When to make your pull request

It is much easier to merge in small changes to the code than extensive ones that touch multiple files. If a small incremental change is possible, please issue a request for that change rather than saving everything up. 

It is a good idea to commit frequently and with commit messages that refer to the changes that have been make and why they have been made. The commit message is usually browsed without seeing the actual changes themselves and which sections of the code have been altered so a helpful message provides relevant details. 

#### Tests

We use `pytest` for the unit testing framework in stripy. In the source directory this means running:

```bash
python setup.py test
```

The existing tests should be passing before you start coding (help us out with an issue if that is not the case !) and when you have finished. Any new functionality should also have tests that we can use to verify the code. It is important that you make it clear if the original tests have had to change to accomodate new code / functionality.

<a name="footnote1">*</a> _Small changes are our favorites as they are much easier to quickly merge into the existing code. Proofreading typos, bug fixes, ... the little stuff is especially welcome._
