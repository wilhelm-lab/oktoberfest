<template>
  <v-btn small color="primary" @click.native="submitTask" v-bind:disabled="submited">
    submit <v-icon>done</v-icon>
  </v-btn>
</template>

<script>
import axios from 'axios'

export default {
  props: ['taskid'],
  data () {
    return {
      submited: false,
      url: process.env.VUE_APP_API_URL + '/api/v1/submitJob',
    }
  },
  methods: {
    showTask(){
      let self = this;
      self.submited = true;
      setTimeout(function() {
        self.$store.commit('reset');
        self.$router.push('/task/' + self.$store.state.task.taskId);
      }, 2000); 
      
    },
    submitTask(){

      let submitObject = this.$store.getters.submitObject;
      let self = this;
      axios.post(
        this.url,
        submitObject,
        { params: { hashInput: this.$store.state.task.hashId } }
      ).then(
        // success
        self.showTask()
      ).catch(function () {
        // error
      });
    }
  }

}
</script>


