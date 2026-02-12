function filterdata(datamat,minMg, minPPi)
    validrows = datamat[:,3].>=minMg .&& datamat[:,4].>=minPPi
    return datamat[validrows,:]
end

function getdatamatrix(filename;minMg = 0, minPPi = 0)
    df = CSV.read(filename,DataFrame)
    filterdata(Matrix(df),minMg, minPPi)
end