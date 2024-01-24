<template>

  <v-container fluid>
  <v-card class='elevation-10 pa-3'>
  <v-card-title class='mt-1 pt-0'>
    <div class="subheading">CSV Format</div>
  </v-card-title>
  <img src="@/assets/csv_example.png" />
  <v-card-text>
    Please provide all three columns below and use <code>,</code> as a separator.
    <ul>
      <li> <code>modified_sequence</code> 
          Use upper case letters in the  column and indicate oxidized Methionine with "M(ox)". 
          Sequence length is restricted to the range of 7 to 30.
          Each C is treated as Cysteine with carbamidomethylation. 
          Prosit does not support U or O as amino acids.
      </li>
      <li> <code>collision_energy</code> 
          Use integer values from 10 and 50.
      </li>
      <li> <code>precursor_charge</code> 
          Use integer values from 1 to 6.
      </li>
      <li> <code>fragmentation</code>
          Either HCD or CID, Use upper case letters.*
      </li>
    </ul>
    <br>
    <b>*Only for TMT model</b>
  </v-card-text>


  </v-card>
  </v-container>
      
   
</template>

<script>
export default {
}
</script>

<style>
  img {
    width: 100%;
    max-width: 800px;
    display: inline-block;
  }
</style>


