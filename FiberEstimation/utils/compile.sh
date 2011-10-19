make -f MakeAverDTI
mv AverDTI ../../bin
rm *.o
make -f MakedtiEigD  
mv dtiEigD ../../bin
rm *.o
make -f MakedtiPD
mv dtiPD ../../bin
rm *.o
make -f MakePDcolorMap
mv PDcolorMap ../../bin
rm *.o
make -f MakeCreateInterleavedDTI
mv createInterleavedDTI ../../bin
rm *.o
make -f MakedtiFA
mv edtiFA ../../bin
rm *.o
make -f MakedtiTrace
mv dtiTrace ../../bin
rm *.o
make -f MakeSmoothDTI
mv smoothDTI ../../bin
rm *.o
