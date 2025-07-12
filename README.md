# Analysis of model predictive control in numbered-up modular facilities

> Code for the following published article:  
> **Yi Dai, Samuel Fay, Andrew Allman**, *“Analysis of model predictive control in numbered-up modular facilities”*, Digital Chemical Engineering, 2023r.  
> DOI: [[10.1016/j.dche.2023.100088]](https://doi.org/10.1016/j.dche.2023.100088)

## 📖 Paper Overview

> Modular production units are a key enabling technology for distributed manufacturing, an important emerging paradigm which helps promote sustainability, resource independence, and robustness in critical chemical
production and distribution infrastructure. In modular systems, desired material throughput is achieved by ‘‘numbering up’’ rather than ‘‘scaling up’’. Previous studies have quantified the benefits of introducing modular
units when making decisions at a supply chain or process design level. However, modular units provide additional potential benefits at the level of operation and control which heretofore has not been as thoroughly
studied. In this paper, we analyze the structure of the optimal modular system control problem and identify three qualities unique to problems with numbered-up units: the ability to control outputs locally or globally,
the presence of symmetry in the underlying problem, and the possibility of reconfiguring the connectivity of modular units. We provide analysis and insights on ways to mitigate challenges and exploit beneficial properties
resulting from the above qualities through a benchmark modular reactor case study.

## ⚙️ Environment & Dependencies

1. **Julia**  
   - Tested with Julia 1.8.1  
2. **Important packages**  
   Ipopt 3.14
   BARON 20.10
   JuMP 0.22.1
