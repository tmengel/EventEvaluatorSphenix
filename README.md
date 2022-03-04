# sPHENIX Event Evaluation
Build the current version via:

cd src
mkdir build
cd build && ./autogen --prefix=$(INSTALLDIR)
make && make install

Run the evaluator on sPHENIX macro:

cd macro
root
.x Fun4All_TannerEval_sPHENIX.C(1)

Designed to run in sphenix environment. Currently reworking to run with sphenixTreeAnalysis
