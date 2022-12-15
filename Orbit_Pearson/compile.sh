# Compile orbit
echo "Compiling Orbit with Pearson Noise..."

g++ -o ct.exe CT.cpp ../house.a -I.. -I./../eigen  -g -Wall -pedantic -O3

g++ -o st.exe statct.cpp -I.. -I./../eigen -g -Wall -pedantic -O3

