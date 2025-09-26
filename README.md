## IIT-Tree: An efficient index to support interval-based query on large temporal graphs
This repository provides the implementation of algorithms for **Interval-based Durability and Existence (IDE) queries** on temporal graphs.  
It includes our proposed methods (`IIT`, `IIT-R`) as well as baselines (`BaseLine`, `DeltaGraph`, `VILA`, `LLAMA`, `DELTA`).

## Code Structure

The main implementations are located in the `query` directory, with each algorithm in its own files:

- `BaseLine.hpp` – Baseline approach
- `IIT.hpp` – Interval Index Tree
- `IIT-R.hpp` – Roaring-compressed Interval Index Tree
- `DeltaGraph.h/cpp` – Delta-based temporal index
- `VILA.h/cpp` – Vectorized Interval Lifespan Algorithm
- `LLAMA.hpp` – Snapshot-based storage with page management
- `DELTA.h/cpp` – Variation-tolerant delta-based structure

---

## Compile the codes
When you already download the codes, run the following commands to compile our codes.
```
mkdir build
cd build
cmake ..
make
```
After running the codes, there will be executable files called `IDE` in `build/query` directory.

---

## Datasets

- **Default datasets** (e.g., `mo3D`, `mo7D`, `mo1M`, `mo2M`, `mo4M`) are included in the package.
- **Additional datasets** can be downloaded from [SNAP](https://snap.stanford.edu/data/index.html).
- Conversion scripts:
    - Use `sim.ipynb` to transform datasets into temporal graph snapshots.

---

## Run the procedure
### Example Commands
- **Run existence query on `IIT`:**
```bash
./query/IDE -d mo1M -m IIT -f 1 -a 2 -b 64
```
- **Run temporal k-core query on `IIT`:**
```bash
./query/IDE -d mo1M -m IIT -f 1 -a 2 -b 64 -q KCORE
```

### Command Line Options

| Option | Required | Description |
|--------|----------|-------------|
| `-d`   | Yes      | Dataset name (e.g., `mo1M`) |
| `-m`   | Yes      | Method to use: `BaseLine`, `IIT`, `IIT-R`, `DeltaGraph`, `VILA`, `LLAMA`, `DELTA` |
| `-f`   | Yes      | Query type: `0` = Durability, `1` = Existence |
| `-a`   | Yes      | Parameter $\alpha$ |
| `-b`   | Yes      | Parameter $\lambda$ |
| `-u`   | No       | Number of update snapshots (**only for IIT**) |
| `-q`   | No       | Specific query type: `KCORE` or `CECI` (**only for IIT**) |

---

## External Dependencies

- **CECI** (used in `IIT`): available at [In-Memory Subgraph Matching: An In-depth Study](https://github.com/RapidsAtHKUST/SubgraphMatching)
- **Peel** (used in `IIT`): available at [kCliqueListing](https://github.com/Gawssin/kCliqueListing)

---