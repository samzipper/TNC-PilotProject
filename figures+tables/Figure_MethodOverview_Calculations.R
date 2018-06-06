require(streamDepletr)

# depletion apportionment
apportion.inv.dist(reach=c(1,2,3,4), dist=c(200, 350, 500, 350), w=2)

# analytical
round(glover(t=365, d=c(200, 350, 500, 350), S=0.1, Tr=1*50), 3)

round(hunt(t=365, d=c(200, 350, 500, 350), S=0.1, Tr=1*50, 
           lmda = streambed_conductance(w=4, Kriv=0.1, briv=1)), 3)

# analytical * depletion apportionment
round(glover(t=365, d=c(200, 350, 500, 350), S=0.1, Tr=1*50)*apportion.inv.dist(reach=c(1,2,3,4), dist=c(200, 350, 500, 350), w=2)[,2], 3)

round(hunt(t=365, d=c(200, 350, 500, 350), S=0.1, Tr=1*50, 
           lmda = streambed_conductance(w=4, Kriv=0.1, briv=1))*
        apportion.inv.dist(reach=c(1,2,3,4), dist=c(200, 350, 500, 350), w=2)[,2], 3)
