function AMPGO02()
   nlp = Model()

   @variable(nlp, x, start=0.0)

   @NLobjective(
    nlp,
    Min,
    sin(x)+sin((10.0/3.0)*x)
   )

   return nlp
end
