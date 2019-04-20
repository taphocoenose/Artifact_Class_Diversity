// This model estimates parameters for the categorical
// distribution of artifact classes across multiple
// regions, based on artifact count data.

data{
   // Number of geographic zones
   int<lower=1> N_zones;
   // Number of artifact classes
   int<lower=1> N_classes;
   // Matrix of artifact class counts, where each row
   // is a geographic zone, and each column is a count
   // of an artifact class.
   int p_counts[N_zones, N_classes];

}
parameters{
   // A vector of parameters, one parameter per artifact
   // class. One vector per geographic zone.
   vector[N_classes] theta_r[N_zones];
   // A scaling parameter for the distribution of
   // theta_r, one per geographic zone, log scale.
   vector[N_zones] theta_scale;
   // Global scaling paramter.
   real scale_mu;
}
model{
   // Prior for log scale scaling parameters
   // across geographic zones
   theta_scale ~ normal(0, 0.5);
   // Prior for the global scaling parameter
   scale_mu ~ normal(1, 0.5);
   
   // Model loop, one iteration for each geographic zone
   for (i in 1:N_zones) {
      // Theta parameters for zone i are normally distributed,
      // with a 0 mean, and an SD specific to zone i.
      theta_r[i] ~ normal(0, exp(theta_scale[i] + scale_mu));
      // Point counts for zone i follow a multinomial
      // distribution. This distribution is defined by
      // a vector of theta parameters calculated from the
      // softmax of zone i's theta_r vector.
      p_counts[i, ] ~ multinomial(softmax(theta_r[i]));
   }
}
generated quantities{
   // Vector of theta parameters, on categorical scale, 
   // for each egographic region.
   simplex[N_classes] theta[N_zones];
   
   // Loop through geographic zones, generating the
   // categorical vector of thetas.
   for (i in 1:N_zones) {
      theta[i] = softmax(theta_r[i]);
   }
   
}
