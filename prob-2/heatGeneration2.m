function heatGen = heatGeneration2(tempIntBd, tempExtBd)

    s = 2;
    hGen = [1; 2; 3; 4; 5];
    avTempA = zeros(5,1);

    for i = 1:length(hGen)
      [numElemWithHG, avTempA(i), globalAvTemp] = poisson(tempIntBd, tempExtBd, hGen(i), 0);
    end

    %plot(avTempA, hGen, 'o-r');
    pHGen = polyfit(avTempA,hGen, 1);
    heatGen = polyval(pHGen, s*tempIntBd);

end
