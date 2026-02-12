# Helpers for loading experimental solubility/precipitation data from CSV.

"""Load a CSV data file as a matrix, optionally filtering rows by minimum Mg and PPi values."""
function getdatamatrix(filename; minMg = 0, minPPi = 0)
    df = CSV.read(filename, DataFrame)
    filterdata(Matrix(df), minMg, minPPi)
end

function filterdata(datamat, minMg, minPPi)
    validrows = datamat[:, 3] .>= minMg .&& datamat[:, 4] .>= minPPi
    return datamat[validrows, :]
end
