# Nonparametric Bayesian multi-armed bandits for single cell experiment design
### Federico Camerlenghi* Bianca Dumitrascu*, Federico Ferrari*, Barbara E. Engelhardt and Stefano Favaro

The problem of maximizing cell type discovery under budget constraints is a fundamental challenge in the collection and the analysis of single-cell RNA-sequencing (scRNA-seq) data. In this paper, we introduce a simple, computationally efficient, and scalable Bayesian nonparametric sequential approach to optimize the budget allocation when designing a large scale collection of scRNA-seq data for the purpose of, but not limited to, creating cell atlases. Our approach relies on i) a hierarchical Pitman-Yor prior that recapitulates biological assumptions regarding cellular differentiation, and ii) a Thompson sampling multi-armed bandit strategy that balances exploitation and exploration to prioritize experiments across a sequence of trials. Posterior inference is performed through a sequential Monte Carlo approach, which allows us to fully exploit the sequential nature of our species sampling problem. We empirically show that our approach outperforms state-of-the-art methods and achieves near-Oracle performance on simulated and real data alike.  HPY-TS code is available at https://github.com/fedfer/HPYsinglecell.


ArXiv: https://arxiv.org/abs/1910.05355
