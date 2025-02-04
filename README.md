# DTALite
DTALite is an open-source, cross-platform, lightweight, and fast Python path engine for networks encoded in [GMNS](https://github.com/zephyr-data-specs/GMNS).


## Quick Start

1. **Tutorial** written in Jupyter notebook with step-by-step demonstration (coming soon).
2. **Documentation** on Installation, Use Cases, Public API, and more (coming soon).


## Installation
DTALite has been published on [PyPI](https://pypi.org/project/DTALite/0.0.1.post1/), and can be installed using
```
$ pip install DTALite
```

### Dependency
The Python modules are written in **Python 3.x**, which is the minimum requirement to explore the most of DTALite.


## Testbed illustration
Users can find the test datasets and code in [test](https://github.com/itsfangtang/DTALite_release/tree/main/test) folder.

**Inputs**: node.csv, link.csv, demand.csv
**Outputs**: link_performance.csv,  od_performance.csv

**The Python code for testing**:
```
import DTALite as dta

dta.dtalite()
```

## How to Cite

Tang, F., Zheng, H., and Zhou, X. (2025, Jan 29). *DTALite*. Retrieved from https://github.com/itsfangtang/DTALite_release


## Please Contribute

Any contributions are welcomed including advise new applications of DTALite, enhance documentation and [docstrings](https://docs.python-guide.org/writing/documentation/#writing-docstrings) in the source code, refactor and/or optimize the source code, report and/or resolve potential issues/bugs, suggest and/or add new functionalities, etc.

Path4GMNS has a very simple workflow setup, i.e., **master for release (on both GitHub and PyPI)** and **dev for development**. If you would like to work directly on the source code (and probably the documentation), please make sure that **the destination branch of your pull request is dev**, i.e., all potential changes/updates shall go to the dev branch before merging into master for release.


## References
Lu, C. C., Mahmassani, H. S., Zhou, X. (2009). [Equivalent gap function-based reformulation and solution algorithm for the dynamic user equilibrium problem](https://www.sciencedirect.com/science/article/abs/pii/S0191261508000829). Transportation Research Part B: Methodological, 43, 345-364.

Jayakrishnan, R., Tsai, W. K., Prashker, J. N., Rajadyaksha, S. (1994). [A Faster Path-Based Algorithm for Traffic Assignment](https://escholarship.org/uc/item/2hf4541x) (Working Paper UCTC No. 191). The University of California Transportation Center.

Bertsekas, D., Gafni, E. (1983). [Projected Newton methods and optimization of multicommodity flows](https://web.mit.edu/dimitrib/www/Gafni_Newton.pdf). IEEE Transactions on Automatic Control, 28(12), 1090–1096.