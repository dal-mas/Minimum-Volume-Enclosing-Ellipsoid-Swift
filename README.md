# Minimum Volume Enclosing Ellipsoid (MVEE)

This code finds the MVEE for a given set of points. It is based on [this](http://www.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid) MatLab algorithm, which is based on the [Khachiyan Algorithm](https://en.wikipedia.org/wiki/Ellipsoid_method).

It uses the [Swix](http://scottsievert.com/swix/) Matrix Library to make the needed matrix operations, so if you are going to use this code, you need to download and install Swix in your project. 