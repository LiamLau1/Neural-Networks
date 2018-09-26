#Neural Networks in C++
Simple Neural Network written in C++ with a test-net which saves the state of the trained neural net

# first_train_net
First neural net courtesy to David Miller, http://millermattson.com/dave
See the associated video for instructions: http://vimeo.com/19569529

#test_net
Takes the trained weights from first-train-net and uses it

#Usage of the simple net
Make sure to change the topology to match the number of inputs in make-training-samples.
Also remember to print out extra inputs
```
g++ -Wall makeTrainingSamples.cpp -o makeTrainingSamples
./makTrainingSamples > /tmp/trainingData.txt
g++ -Wall first_train_net.cpp -o first_train_net
./first_train_net > tmpfile
g++ -Wall predtraining.cpp -o predtraining
./predtraining > /tmp/PredictData.txt
g++ -Wall test_net.cpp -o test_net
./test_net > tmpfile2
```



