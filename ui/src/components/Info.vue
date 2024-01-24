<template>

  <v-container grid-list-md>
    <v-layout row wrap>
      <v-card flat>
        <v-card-text>
          <slot />
        </v-card-text>
      </v-card>
    </v-layout>
  </v-container>
      
   
</template>

<script>
export default {
  props: ['text'],
}
</script>


