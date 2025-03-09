# DTALite
DTALite is an open-source, cross-platform, lightweight, and fast Python path engine for networks encoded in [GMNS](https://github.com/zephyr-data-specs/GMNS).


## Quick Start

1. **[Tutorial](https://github.com/itsfangtang/DTALite_release/blob/main/dtalite_validate.ipynb)** written in Jupyter notebook with step-by-step demonstration.
2. **[Documentation](https://github.com/itsfangtang/DTALite_release/wiki/DTALite-Inputs-and-Outputs)** on inputs preparation and outputs explanations.


## Installation
DTALite has been published on [PyPI](https://pypi.org/project/DTALite/0.0.1.post1/), and can be installed using
```
$ pip install DTALite
```

### Dependency
The Python modules are written in **Python 3.x**, which is the minimum requirement to explore the most of DTALite.


## Testbed illustration
Users can find the test datasets and code in [test](https://github.com/itsfangtang/DTALite_release/tree/main/test) folder.

**Inputs**: node.csv, link.csv, demand.csv, settings.csv
**Outputs**: link_performance.csv,  od_performance.csv

**The Python code for checking and sorting node.csv and link.csv**:
```
def check_and_sort_files(node_file='node.csv', link_file='link.csv'):
    """
    Check if node_file is sorted by node_id and link_file is sorted first by
    from_node_id and then by to_node_id. If not, sort the files and save them.

    Parameters:
    - node_file: str, path to the node CSV file.
    - link_file: str, path to the link CSV file.
    """

    # Check if node.csv is sorted by node_id
    node_df = pd.read_csv(node_file)
    if not node_df['node_id'].is_monotonic_increasing:
        print(f"{node_file} is not sorted by node_id. Sorting now.")
        node_df.sort_values('node_id', inplace=True)
        node_df.to_csv(node_file, index=False)
    else:
        print(f"{node_file} is already sorted by node_id.")

    # Check if link.csv is sorted by from_node_id and to_node_id (first sort by from_node_id, then to_node_id)
    link_df = pd.read_csv(link_file)
    sorted_link_df = link_df.sort_values(by=['from_node_id', 'to_node_id']).reset_index(drop=True)
    current_link_order = link_df[['from_node_id', 'to_node_id']].reset_index(drop=True)

    if not current_link_order.equals(sorted_link_df[['from_node_id', 'to_node_id']]):
        print(f"{link_file} is not sorted by from_node_id and to_node_id. Sorting now.")
        sorted_link_df.to_csv(link_file, index=False)
    else:
        print(f"{link_file} is already sorted by from_node_id and to_node_id.")
```

**The Python code for testing**:
```
import DTALite as dta

# Use the function to check and sort the files, then run the assignment.
check_and_sort_files()
dta.assignment()

# Smulation if needed
## dta.simulation()
```


## How to Cite

Tang, F., Zheng, H., and Zhou, X. (2025, Jan 29). *DTALite*. Retrieved from https://github.com/itsfangtang/DTALite_release


## Please Contribute

Welcome to join the [DTALite Google Group](https://groups.google.com/g/dtalite)! Any contributions are welcomed including advise new applications of DTALite, enhance documentation and [docstrings](https://docs.python-guide.org/writing/documentation/#writing-docstrings) in the source code, refactor and/or optimize the source code, report and/or resolve potential issues/bugs, suggest and/or add new functionalities, etc.

DTALite has a very simple workflow setup, i.e., **main branch for release (on both GitHub and PyPI)** and **dev for development**. If you would like to work directly on the source code (and probably the documentation), please make sure that **the destination branch of your pull request is dev**, i.e., all potential changes/updates shall go to the dev branch before merging into master for release.


## References

**1. White Paper and Application**

Zhou, X., & Taylor, J. (2014). [DTALite: A queue-based mesoscopic traffic
simulator for fast model evaluation and
calibration.](https://www.tandfonline.com/doi/full/10.1080/23311916.2014.961345)
Cogent Engineering 1.1, 961345.

Marshall, N.L. (2018). [Forecasting the impossible: The status quo of estimating
traffic flows with static traffic assignment and the future of dynamic traffic
assignment.](https://www.sciencedirect.com/science/article/pii/S2210539517301232)
Research in Transportation Business & Management 29, 85-92.

**2. [NeXTA/DTALite Workshop Webinar](https://www.youtube.com/channel/UCUHlqojCQ4f7VvqroUhbaFA)**

**3. [Mini-Lesson on the Algorithmic Details](https://youtu.be/rorZAhNNOf0)**

What is the best way to learn dynamic traffic simulation and network assignment
for a beginner? Do you want to integrate a powerful traffic simulator in your deep
learning framework? We would like to offer a collaborative learning experience
through 500 lines of Python codes and real-life data sets. This is part of our
mini-lessons through teaching dialog.

**4. Parallel Computing Algorithms**

Qu, Y., & Zhou, X. (2017). Large-scale dynamic transportation network simulation:
A space-time-event parallel computing approach. Transportation Research Part C:
Emerging technologies, 75, 1-16.

**5. OD Demand Estimation**

Lu, C. C., Zhou, X., & Zhang, K. (2013). Dynamic origin–destination demand flow
estimation under congested traffic conditions.Transportation Research Part C:
Emerging Technologies, 34, 16-37.

**6. Simplified Emission Estimation Model**

Zhou, X., Tanvir, S., Lei, H., Taylor, J., Liu, B., Rouphail, N. M., & Frey, H. C. (2015). Integrating a simplified emission estimation model and mesoscopic dynamic traffic simulator to efficiently evaluate emission impacts of traffic management strategies.Transportation Research Part D: Transport and Environment, 37, 123-136.

**7. Eco-system Optimal Time-dependent Flow Assignment**

Lu, C. C., Liu, J., Qu, Y., Peeta, S., Rouphail, N. M., & Zhou, X. (2016). Eco-system optimal time-dependent flow assignment in a congested network. Transportation Research Part B: Methodological, 94, 217-239.

**8. Transportation-induced Population Exposure Assessment**

Vallamsundar, S., Lin, J., Konduri, K., Zhou, X., & Pendyala, R. M. (2016). A
comprehensive modeling framework for transportation-induced population exposure
assessment. Transportation Research Part D: Transport and Environment, 46, 94-113.

**9. Integrated ABM and DTA**

Xiong, C., Shahabi, M., Zhao, J., Yin, Y., Zhou, X., & Zhang, L. (2020). An
integrated and personalized traveler information and incentive scheme for energy
efficient mobility systems. Transportation Research Part C: Emerging Technologies,
113, 57-73.

**10. State-wide Transportation Modeling**

Zhang, L. (2017). Maryland SHRP2 C10 Implementation Assistance – MITAMS: Maryland
Integrated Analysis Modeling System, Maryland State Highway Administration.

**11. Workzone Applications**

Schroeder, B, et al. (2014). Work zone traffic analysis & impact assessment. FHWA/NC/2012-36.
North Carolina. Dept. of Transportation. Research and Analysis Group.