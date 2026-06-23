# Automated Lower Bounds for Bilinear Complexity over Finite Fields

We present a general, automated framework for proving lower bounds on the
bilinear complexity (tensor rank) of multiplication problems over a finite field
$`\mathbb{F}_q`$. The framework is parameterized only by the multiplication tensor and by a
group of rank-preserving symmetries acting on one argument: it classifies
the orbits of constraint subspaces under that group, runs a dynamic program over
the orbits combining four lower-bound techniques, and emits a proof certificate
that a verifier rechecks, typically faster than the search.

Instantiating the framework for matrix multiplication, we improve the lower bounds
for four small formats over $`\mathbb{F}_2`$, most notably showing that
the bilinear complexity of multiplying two $`3 \times 3`$ matrices over $`\mathbb{F}_2`$ is at
least $`20`$, raising the bound of $`19`$ that had stood since Bläser (2003).
Instantiating it for polynomial multiplication — full products,
cyclic convolution, and the truncated (modulo $`x^N`$) and negacyclic (modulo
$`x^N+1`$) products — we obtain eighteen new lower bounds over $`\mathbb{F}_2`$ and $`\mathbb{F}_3`$.
Every bound is backed by a machine-checkable certificate.

For full details, definitions, and proofs, see the [paper](https://arxiv.org/abs/2603.07280).

## Results

### Matrix multiplication

The framework improves the best known lower bound for four small formats over
$`\mathbb{F}_2`$:

$$\mathbf{R}_{\mathbb{F}_2}(\langle 2,3,4\rangle)\ge 19,\quad
  \mathbf{R}_{\mathbb{F}_2}(\langle 3,3,3\rangle)\ge 20,\quad
  \mathbf{R}_{\mathbb{F}_2}(\langle 3,3,4\rangle)\ge 25,\quad
  \mathbf{R}_{\mathbb{F}_2}(\langle 3,4,4\rangle)\ge 29.$$

### Polynomial multiplication

The same framework yields eighteen new lower bounds over $`\mathbb{F}_2`$ and
$`\mathbb{F}_3`$, across the full product and the three quotient families
$`\mathbb{F}_q[x]/(x^N-\gamma)`$. Write $`\mathbf{R}_{\mathbb{F}_q}(\mathrm{P}_N)`$
for the bilinear complexity of the full product of two degree-$`(N-1)`$
polynomials, and $`\mathbf{R}_{\mathbb{F}_q}(\mathrm{C}_N)`$,
$`\mathbf{R}_{\mathbb{F}_q}(\mathrm{T}_N)`$, and
$`\mathbf{R}_{\mathbb{F}_q}(\mathrm{C}^{-}_N)`$ for, respectively, cyclic
convolution (modulo $`x^N-1`$), the truncated product (modulo $`x^N`$), and the
negacyclic product (modulo $`x^N+1`$), each of length $`N`$. Over $`\mathbb{F}_2`$ the
negacyclic product coincides with the cyclic one, since $`x^N+1=x^N-1`$.

Full product:

$$\mathbf{R}_{\mathbb{F}_2}(\mathrm{P}_6)\ge 16,\quad
  \mathbf{R}_{\mathbb{F}_2}(\mathrm{P}_7)\ge 19,\quad
  \mathbf{R}_{\mathbb{F}_2}(\mathrm{P}_8)\ge 21,\quad
  \mathbf{R}_{\mathbb{F}_3}(\mathrm{P}_6)\ge 14.$$

Cyclic convolution:

$$\mathbf{R}_{\mathbb{F}_2}(\mathrm{C}_7)\ge 13,\quad
  \mathbf{R}_{\mathbb{F}_2}(\mathrm{C}_8)\ge 19,\quad
  \mathbf{R}_{\mathbb{F}_2}(\mathrm{C}_{10})\ge 22,\quad
  \mathbf{R}_{\mathbb{F}_3}(\mathrm{C}_9)\ge 19.$$

Truncated product:

$$\mathbf{R}_{\mathbb{F}_2}(\mathrm{T}_6)\ge 13,\quad
  \mathbf{R}_{\mathbb{F}_2}(\mathrm{T}_7)\ge 16,\quad
  \mathbf{R}_{\mathbb{F}_2}(\mathrm{T}_8)\ge 19,$$

$$\mathbf{R}_{\mathbb{F}_2}(\mathrm{T}_9)\ge 21,\quad
  \mathbf{R}_{\mathbb{F}_2}(\mathrm{T}_{10})\ge 24,\quad
  \mathbf{R}_{\mathbb{F}_2}(\mathrm{T}_{11})\ge 27,$$

$$\mathbf{R}_{\mathbb{F}_3}(\mathrm{T}_8)\ge 17,\quad
  \mathbf{R}_{\mathbb{F}_3}(\mathrm{T}_9)\ge 19,\quad
  \mathbf{R}_{\mathbb{F}_3}(\mathrm{T}_{10})\ge 21.$$

Negacyclic product:

$$\mathbf{R}_{\mathbb{F}_3}(\mathrm{C}^{-}_9)\ge 19.$$

## Build

This Git repository requires [LFS](https://git-lfs.com/). Install it before `git clone`.

We use Bazel (via Bazelisk) to build the project. Install Bazelisk from `https://github.com/bazelbuild/bazelisk`.

This project requires C++20 support. Configure your C++ toolchain in `.bazelrc`.

```bash
# release build
bazel build --config=opt //...

# run the unit tests
bazel test --config=opt //...
```

## Running the pipeline (`run.py`)

The subcommand `verify` verifies a single certificate. For example,
```bash
python3 run.py verify certs/matrix/cert_matrix_q02_n333.pb.txt
```
verifies the tensor rank lower bounds for the 3x3x3 matrix multiplication problem over $`\mathbb{F}_2`$.

Other subcommands of `run.py` tell how `certs/` are generated.

## Code layout

The directory split mirrors the verifier's trust boundary — the set of code a
reviewer must trust to believe a reported rank **lower** bound. See
[`verifier/README.md`](verifier/README.md).

- **Trust base (audit this):** `verifier/` (the checker), `core/` (field
  arithmetic, tensors, constraints, the certificate proto and proof codecs), and
  `problems/` (the tensor and symmetry-group definitions, plus the build-time
  problem selection).
- **Not trusted:** `search/` (orbit enumeration, the meet-in-the-middle orbit
  map, and the lower-bound proof *generators*) and `upper_bound/` (the flip-graph
  rank upper bound). A bug here can only make the verifier *reject*, never wrongly
  accept.

Four problem families live under `problems/`: full and cyclic polynomial
multiplication, extension-field multiplication, and ⟨N0,N1,N2⟩ matrix
multiplication.

## Citation

```bibtex
@misc{wang2026automatedlowerboundstensor,
      title={Automated Lower Bounds for Bilinear Complexity over Finite Fields}, 
      author={Chengu Wang},
      year={2026},
      eprint={2603.07280},
      archivePrefix={arXiv},
      primaryClass={cs.CC},
      url={https://arxiv.org/abs/2603.07280}, 
}
```
