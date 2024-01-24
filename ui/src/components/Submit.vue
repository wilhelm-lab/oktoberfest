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
      url: '/prosit/api/submit.xsjs',
    }
  },
  methods: {
    showTask(){
      let self = this;
      setTimeout(function() {
        let taskid = self.taskid;
        self.$store.commit('reset');
        self.$router.push('/prosit/task/' + taskid);
      }, 2000); 
      
    },
    submitTask(){

      let submitObject = this.$store.getters.submitObject;
      let self = this;
      this.submited = true;
      submitObject.id = this.taskid;
      axios.post(
        this.url,
        submitObject,
        { params: { datasetId: this.taskid } }
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


