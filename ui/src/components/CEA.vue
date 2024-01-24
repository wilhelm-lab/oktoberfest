<template>

<v-stepper v-model="step" vertical>

  <Info> 
      This task estimates the optimal collision energy (CE) based on a given search result.
      You need to upload a RAW file as well as the MaxQuant's msms.txt for calibration. <br />
      Prosit will: 
      <ol>
        <li>Select a random subset of high-scoring PSMs</li>
        <li>Predict those in for each CE from 18 to 39.</li>
        <li>Calculate which CE achieves highest correlations with the experimental spectra</li>
      </ol>
  Please note: Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported.
  Each C is treated as Cysteine with carbamidomethylation (fixed modification in MaxQuant).<br />
  Also note: Antivirus software may cancel large uploads - turn it off if you experience upload resets.
  </Info>


  <v-stepper-step :complete="step > 1" step="1">
    Upload Files
    <small>msms.txt and RAW file.</small>
  </v-stepper-step>
  <v-stepper-content step="1">
    
    <Upload filetype="msms.txt" filesuffix=".txt" :taskid="taskid" hinttext="MaxQuant's msms.txt from a finished search. Note, amino acid U or O are not supported"></Upload>
    <Upload filetype="RAW" filesuffix=".RAW,.raw"  :taskid="taskid" hinttext="RAW file that was searched (restricted to Thermo Fisher HCD Orbitrap). File size is limited to 2GB."></Upload>

    <v-card flat>   
    <v-card-actions>
    <v-spacer></v-spacer>
    <v-btn small color="primary" @click="step = 2" :disabled="!doneUploads">next<v-icon>chevron_right</v-icon></v-btn>
    </v-card-actions>
    </v-card>
  </v-stepper-content>

   <v-stepper-step :complete="step > 2" step="2">
    Model
    <small>Select intensity and iRT model for prediction</small>
  </v-stepper-step>  
  <v-stepper-content step="2">
  <v-card flat>
    <div v-if="modelIntensityList.length>0">
      Intensity prediction model
      <v-radio-group v-model="modelIntensityName" column>
        <v-radio color='primary'
                 v-for="f in modelIntensityList"
                 :key="f.name"
                 :label="f.name"
                 :value="f.name"
        ></v-radio>
      </v-radio-group>
    </div>

    <div v-if="modelIRTList.length>0">
      iRT prediction model
      <v-radio-group v-model="modelIRTName" column>
        <v-radio color='primary'
                 v-for="f in modelIRTList"
                 :key="f.name"
                 :label="f.name"
                 :value="f.name"
        ></v-radio>
      </v-radio-group></div>

  <v-card-actions>
  <v-spacer></v-spacer>
  <v-btn small color="accent" @click="step = 1"><v-icon>chevron_left</v-icon>back</v-btn>
  <v-btn small color="primary" @click="toStepThree" :disabled="!((modelIRTName||modelIRTList.length==0)&&(modelIntensityName||modelIntensityList.length==0))">next<v-icon>chevron_right</v-icon></v-btn>
  </v-card-actions>
  </v-card> 
  </v-stepper-content>

  <v-stepper-step :complete="step > 3" step="3">
    Isobaric Label
    <small>If TMT model is selected</small>
  </v-stepper-step>
  <v-stepper-content step="3">
    <v-card flat>
      <v-radio-group v-model="tag" column>
        <v-radio color='primary'
                 label="TMT6/10/11-plex"
                 value="tmt"
        ></v-radio>
        <v-radio color='primary'
                 label="iTRAQ4-plex"
                 value="itraq4"
        ></v-radio>
        <v-radio color='primary'
                 label="iTRAQ8-plex"
                 value="itraq8"
        ></v-radio>
        <v-radio color='primary'
                 label="TMT16/18-plex (TMTPro)"
                 value="tmtpro"
        ></v-radio>
      </v-radio-group>

      <v-card-actions>
        <v-spacer></v-spacer>
        <v-btn small color="accent" @click="step = 2"><v-icon>chevron_left</v-icon>back</v-btn>
        <v-btn small color="primary" @click="toStepFour" :disabled="!((tag))">next<v-icon>chevron_right</v-icon></v-btn>
      </v-card-actions>
    </v-card>
  </v-stepper-content>

  <v-stepper-step :complete="step > 4" step="4">
    Task ID
    <small>Check if everything is correct and submit the task</small>
  </v-stepper-step>  
  <v-stepper-content step="4">
  <v-card flat>
  <v-card-actions>
  <v-spacer></v-spacer>
  <v-btn small color="accent" @click="step = 3"><v-icon>chevron_left</v-icon>back</v-btn>
  <Submit :taskid="taskid"/>
  </v-card-actions>
  </v-card> 
  </v-stepper-content>
        
</v-stepper>


</template>

<script>
import axios from 'axios'
import Upload from '@/components/Upload.vue'
import Submit from '@/components/Submit.vue'
import Info from '@/components/Info.vue'


export default {
  props: ['taskid'],
  components: { Upload, Submit, Info },
  data () {
    return {
      step: 1,
      modelIntensityName: false,
      modelIntensityList: false,
      modelIRTName: false,
      modelIRTList: false,
      tag: false,
    }
  },
  computed: {
    doneUploads: function () {
      let uploads = this.$store.state.task.uploads;
      return uploads['msms.txt'] && uploads.RAW;
    }
  },
  methods: {
    toStepThree: function (){
      this.$store.commit('changeIRTModel', this.modelIRTName);
      this.$store.commit('changeIntensityModel', this.modelIntensityName);
      if(this.modelIRTName.includes('TMT') || this.modelIntensityName.includes('TMT')){
        this.step = 3;
      }else{
        this.step = 4;
      }
    },
    toStepFour: function (){
      this.$store.commit('changeTag', this.tag);
      this.step = 4;
    },
    setModels: function(response){
      let modelIntensityLists = response.data.models.intensities;
      let modelIRTLists = response.data.models.iRT;
      this.modelIRTList = [];
      //modelIRTList = response.data.models.iRT;
      this.modelIntensityList = [];
      for (let i = 0; i < modelIntensityLists.length; i++){
        //if(modelIntensityLists[i].enabled.ce_calibration)
        this.modelIntensityList.push(modelIntensityLists[i]);
      }
      for (let i = 0; i < modelIRTLists.length; i++){
        //if(modelIRTLists[i].enabled.ce_calibration)
        this.modelIRTList.push(modelIRTLists[i]);
      }
      let that = this;
      if(this.modelIntensityList.length>0) {
        this.modelIntensityList.forEach(function (item) {
          if (item.default) that.modelIntensityName = item.name;
        });
      }
      if(this.modelIRTList.length>0) {
        this.modelIRTList.forEach(function (item) {
          if (item.default) that.modelIRTName = item.name;
        });
      }
    }
  },
  mounted() {
    let modelsUrl = '/prosit/api/models.xsjs';
    axios.get(modelsUrl).then(this.setModels);
  }
  
}
</script>

<style scoped></style>