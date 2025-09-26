## Compile the codes
When you already download the codes, run the following commands to compile our codes.
```
mkdir build
cd build
cmake ..
make
```
After running the codes, there will be executable files called `IDE`, which means you have already compiled our codes. The main code of BaseLine, IIT, DeltaGraph and VILA are in `main.cpp`.
## Run the procedure
```
./bit_test -d mo4M -m 1 -f 1 -k 2 -e 64
```
|Syntax|Description|
|---|---|
|-d|Dataset|
|-m|Method(BaseLine 0, ITT 1, ITT-R-2, DeltaGraph 3, VILA 4)|
|-f|Durability 0 Or Existence 1|
|-k|$\alpha$|
|-e|$\lambda$|

The others datasets are in [SNAP](https://snap.stanford.edu/data/index.html), the code to convert dataset to snapshots is called `sun.ipynb`.
 
