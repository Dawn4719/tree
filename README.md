## Compile the codes
When you already download the codes, run the following commands to compile our codes.
```
mkdir build
cd build
cmake ..
make
```
After running the codes, there will be an executable files called `bit_test`, which means you have already compiled our codes.
## Run the procedure
```
./bit_test -d su4M -m 1 -f 1 -k 2 -e 64
```
|Syntax|Description|
|---|---|
|-d|Dataset|
|-m|Method(BaseLine 0, ITT 1, DeltaGraph 2, VILA 3)|
|-f|Durability Or Existence|
|-k|$\alpha$|
|-e|$\lambda$|
 
