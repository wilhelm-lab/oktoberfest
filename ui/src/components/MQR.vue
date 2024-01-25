<template>

<v-stepper v-model="step" vertical>

  <Info> 
      This task rescores an existing MaxQuant search (FDR 100%) using features generated from fragmentation prediction. 
      You need to upload a RAW file as well as the MaxQuant's msms.txt file from a search. <br />
      Prosit will: 
      <ol>
        <li>Calibrate itself against the RAW.</li>
        <li>Predict all sequences in the msms.txt.</li>
        <li>Use the predicted spectra to generate features for percolator.</li>
        <li>Run percolator to rescore the search. </li>
      </ol>
      Please note: You need a MaxQuant search at 100% FDR, otherwise targets may be filtered by MaxQuant's FDR calculation before rescoring. 
      Sequences with amino acid U or O are not supported. Modifications except "M(ox)" are not supported.
      Each C is treated as Cysteine with carbamidomethylation (fixed modification in MaxQuant).<br />
      Also note: Antivirus software may cancel large uploads - turn it off if you experience upload resets.
  </Info>



  <v-stepper-step :complete="step > 1" step="1">
    Search Engine
    <small>select the search engine used for searching your data.</small>
  </v-stepper-step>
  <v-stepper-content step="1">
     <v-card flat>   
      <div v-if="searchEngineList.length>0">
        Intensity prediction model
        <v-radio-group v-model="searchEngine" column>
          <v-radio color='primary'
                  v-for="f in searchEngineList"
                  :label="f"
                  :value="f"
          ></v-radio>
        </v-radio-group>
      </div>
    <v-card-actions>
    <v-spacer></v-spacer>
    <v-btn small color="primary" @click="toStepTwo" :disabled="!((searchEngine))">next<v-icon>chevron_right</v-icon></v-btn>
    </v-card-actions>
    </v-card>
  </v-stepper-content>

  <v-stepper-step :complete="step > 2" step="2">
    Upload Files
    <div v-if="searchEngine='MaxQuant'">
      <small>msms.txt and RAW file.</small>
    </div>
    <div v-if="searchEngine='MSFragger'">
      <small>pepXML files and RAW file.</small>
    </div>
    <div v-if="searchEngine='Sage'">
      <small>psm.tsv and RAW file.</small>
    </div>
  </v-stepper-step>
  <v-stepper-content step="2">
    <div v-if="searchEngine='MaxQuant'">
      <Upload filetype="msms.txt" filesuffix=".txt" :taskid="taskid" hinttext="MaxQuant's msms.txt from a finished search. Note, amino acid U or O are not supported."></Upload>
    </div>
    <div v-if="searchEngine='MSFragger'">
      <MultiUpload filetype="pepXML" filesuffix=".pepXML" :taskid="taskid" hinttext="MSFragger .pepXML files generated during the search"></MultiUpload>
    </div>
    <div v-if="searchEngine='Sage'">
      <Upload filetype="results.sage.tsv" filesuffix=".tsv" :taskid="taskid" hinttext="Sage search results."></Upload>
    </div>
    <MultiUpload filetype="RAW" filesuffix=".RAW,.raw,.zip" v-bind:sizelimit=6000 :taskid="taskid" hinttext="RAW file that was searched."></MultiUpload>

    <v-card flat>   
    <v-card-actions>
    <v-spacer></v-spacer>
    <v-btn small color="primary" @click="step = 3" :disabled="!doneUploads">next<v-icon>chevron_right</v-icon></v-btn>
    </v-card-actions>
    </v-card>
  </v-stepper-content>

  <v-stepper-step :complete="step > 3" step="3">
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
  <v-btn small color="accent" @click="step = 2"><v-icon>chevron_left</v-icon>back</v-btn>
  <v-btn small color="primary" @click="toStepFour" :disabled="!((modelIRTName||modelIRTList.length==0)&&(modelIntensityName||modelIntensityList.length==0))">next<v-icon>chevron_right</v-icon></v-btn>
  </v-card-actions>
  </v-card> 
  </v-stepper-content>


  <v-stepper-step :complete="step > 4" step="4">
    Isobaric Label
    <small>If TMT model is selected</small>
  </v-stepper-step>
  <v-stepper-content step="4">
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
        <v-btn small color="accent" @click="step = 3"><v-icon>chevron_left</v-icon>back</v-btn>
        <v-btn small color="primary" @click="toStepFive" :disabled="!((tag))">next<v-icon>chevron_right</v-icon></v-btn>
      </v-card-actions>
    </v-card>
  </v-stepper-content>

  <v-stepper-step :complete="step > 5" step="5">
    Task ID
    <small>Check if everything is correct and submit the task</small>
  </v-stepper-step>
  <v-stepper-content step="5">
    <v-card flat>
      <v-card-actions>
        <v-spacer></v-spacer>
        <v-btn small color="accent" @click="step = 5"><v-icon>chevron_left</v-icon>back</v-btn>
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
import MultiUpload from '@/components/MultiUpload.vue'


export default {
  props: ['taskid'],
  components: { Upload, Submit, MultiUpload, Info},
  data () {
    return {
      step: 1,
      modelIntensityName: false,
      modelIntensityList: false,
      searchEngineList: false,
      modelIRTName: false,
      modelIRTList: false,
      tag: false,
      searchEngine: false,
    }
  },
  computed: {
    doneUploads: function () {
      let uploads = this.$store.state.task.uploads;
      return uploads['msms.txt'] && uploads.RAW;
    }
  },
  methods: {
    toStepTwo: function (){
      this.$store.commit('changeSearchEngine', this.searchEngine);
      this.step = 2;
    },
    toStepFour: function (){
      this.$store.commit('changeIRTModel', this.modelIRTName);
      this.$store.commit('changeIntensityModel', this.modelIntensityName);
      if(this.modelIRTName.includes('TMT') || this.modelIntensityName.includes('TMT')){
        this.step = 4;
      }else{
        this.tag='tmt';
        this.step = 5;
      }
    },
    toStepFive: function (){
      this.$store.commit('changeTag', this.tag);
      this.step = 5;
    },
    setModels: function(response){
      let modelIntensityLists = response.data.models.intensities;
      let modelIRTLists = response.data.models.iRT;
      this.modelIRTList = [];
      this.modelIntensityList = [];
      this.searchEngineList = ['MaxQuant','MSFragger','Sage'];
      for (let i = 0; i < modelIntensityLists.length; i++){
        //if(modelIntensityLists[i].enabled.rescoring)
        this.modelIntensityList.push(modelIntensityLists[i]);
      }
      for (let i = 0; i < modelIRTLists.length; i++){
        //if(modelIRTLists[i].enabled.rescoring)
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