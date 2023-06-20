clearvars
close all
%clc

tempIntBd = 55.0;
tempExtBd = 85.0;
R1 = 39.9;
R2 = 60.1;
s = 2;
rHeatGen = 124.9;

hGen = [1; 2; 3; 4; 5];
avTempA = zeros(5,1);

for i = 1:length(hGen)
   [numElemWithHG, avTempA(i), globalAvTemp] = poisson(tempIntBd, tempExtBd, hGen(i), 0);
end

plot(avTempA, hGen, 'o-r')
pHGen = polyfit(avTempA,hGen, 1);
heatGenInterp = polyval(pHGen, s*tempIntBd);

% Solutions
fprintf('Thermal problem. Part (c)\n')
fprintf('Temperature at internal boundary, tempIntBd = %e\n', tempIntBd)
fprintf('Temperature at external boundary, tempExtBd = %e\n', tempExtBd)
fprintf('Number of elements with heat generation: %d\n', numElemWithHG)
fprintf('Heat generation at R > %f, using interpolation, f = %e\n',...
    rHeatGen, heatGenInterp)

heatGenLinearSystem = heatGeneration1(tempIntBd, tempExtBd, 0);

fprintf('Heat generation at R > %f, using linear system, f = %e\n',...
    rHeatGen, heatGenLinearSystem)

err = abs(heatGenInterp-heatGenLinearSystem);

fprintf('Difference between the two methods of computing heatGen, ')
fprintf('err = |heatGenInterp - heatGenLinearSystem| = %.4e\n',err)

% Check Solution
[numElemWithHeatGen, avTempAnnulus, globalAvTemp] = ...
    poisson(tempIntBd, tempExtBd, heatGenInterp, 1);

% Solutions
fprintf('Averaged temprature in the annulus %f < R < %f, <T> = %.4e\n',...
    R1, R2, avTempAnnulus)
fprintf('Err:= |<T at R> - %d x intBdTemp| = |<T at R> - %f | = %.4e\n',...
    s, s*tempIntBd, abs(avTempAnnulus-s*tempIntBd))
fprintf('Hint: global Averaged temperature, <T> = %.4e\n',  globalAvTemp)
