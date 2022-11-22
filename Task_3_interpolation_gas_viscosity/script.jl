using Visc_new

CO_2 = Substance(44.009, 3.996, 190.0, 300.0, 1000.0)
CH_4 = Substance(16.043, 3.822, 137.0, 100.0, 600.0)
O_2 =  Substance(31.999, 3.433, 113.0, 100.0, 1000.0)
#println("eta: ", viscos(500, CH_4))
eta = zeros(Int(O_2.T_max))

for T in Int(O_2.T_min):Int(O_2.T_max)
    eta[T] = viscos(T, O_2)
end
print(eta)
