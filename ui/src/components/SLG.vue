<template>

<v-stepper v-model="step" vertical>

  <Info> 
      This task generates a spectral library either by digesting a given FASTA file, 
      or by predicting a list of peptides given in a CSV file. 
      You need to provide a collision energy (CE) for prediction. 
      To estimate an optimal CE for prediction, please use "CE Calibration".
       <br />
      When a FASTA file is provided, Prosit will: 
      <ol>
        <li>Digest the FASTA, for the given parameters (i.e. protease).</li>
        <li>Predict all spectra at the given collision energy.</li>
      </ol>
      When a CSV with peptides is provided, Prosit will directly predict all spectra.<br />
      Please note: Antivirus software may cancel large uploads - turn it off if you experience upload resets.
  </Info>



  <v-stepper-step :complete="step > 1" step="1">
    Settings
    <small>Indicate collision energy, the maximum number of missed cleavages, and number of oxidized methionines per peptide.</small>
  </v-stepper-step>
  <v-stepper-content step="1">   
    <v-card flat>
        How would you like to provide the list of peptides?

        <v-radio-group v-model="fromFasta" >
          <v-radio
            color='primary'
            key="1"
            label="CSV"
            :value="false"
          ></v-radio>
          
          <v-radio disabled
            color='primary'
            key="2"
            label="FASTA (comming soon)"
            :value="true"
          ></v-radio>
          
        </v-radio-group>

        <CSV v-show="!fromFasta"/>
   
        <v-container>
        <v-card flat v-show="fromFasta" class='elevation-10 pa-3'>
          <v-card-title class='mt-1 pt-0'>
            <div class="subheading">FASTA digestion parameters</div>
          </v-card-title>
          <v-card-text>
        <v-slider prepend-icon="blur_linear" v-model="collisionEnergy" step=1 min=9 max=40 thumb-label='always' label="collision energy"></v-slider>
        <v-select prepend-icon="crop" :items="proteases" v-model="protease" label="protease" ></v-select>
        <v-select prepend-icon="tune" :items="missedClevages" v-model="maxMissedClevages" label="missed cleavages (max)"></v-select>
        <v-select prepend-icon="toll" :items="oxidizedMethionine" v-model="maxOxidizedMethionine" label="oxidized methionine (max)"></v-select>
        </v-card-text>
        </v-card>
        </v-container>
        
      <v-card-actions>
      <v-spacer></v-spacer>
      <v-btn small color="primary" @click="toStepTwo">next<v-icon>chevron_right</v-icon></v-btn>
      </v-card-actions>
    </v-card>
  </v-stepper-content>


  <v-stepper-step :complete="step > 2" step="2">
    Upload Files
    <small>Fasta or CSV with list of peptides</small>
  </v-stepper-step>
  <v-stepper-content step="2">

    <CSV v-show="!fromFasta"/>
    <Upload :hidden="!fromFasta" filetype="FASTA" filesuffix=".FASTA,.fasta" :taskid="taskid" hinttext="FASTA to be digested to a peptide list."  />
    <Upload :hidden="fromFasta" filetype='peptides.csv' filesuffix=".csv" :taskid="taskid" hinttext="containing peptide sequence, collision energy and precursor charge." />

    <v-card flat>   
    <v-card-actions>
    <v-spacer></v-spacer>
    <v-btn small color="accent" @click="step = 1"><v-icon>chevron_left</v-icon>back</v-btn>
    <v-btn small color="primary" @click="toStepThree" :disabled="!doneUploads">next<v-icon>chevron_right</v-icon></v-btn>
    </v-card-actions>
    </v-card> 
  </v-stepper-content>

   <v-stepper-step :complete="step > 3" step="3">
    Model
    <small>Select intensity and iRT model for prediction</small>
  </v-stepper-step>  
  <v-stepper-content step="3">
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
    <div :hidden="!fromFasta">
      <v-tooltip top>
        <v-chip label outline color="primary" slot="activator">
          <v-icon left>blur_linear</v-icon>{{collisionEnergy}}
        </v-chip>
        <span>collision energy</span>
      </v-tooltip>
      <v-tooltip top>
        <v-chip label outline color="primary" slot="activator">
          <v-icon left>crop</v-icon>{{protease}}
        </v-chip> 
        <span>protease</span>
      </v-tooltip>
      <v-tooltip top>
        <v-chip label outline color="primary" slot="activator">
            <v-icon left>tune</v-icon>{{maxMissedClevages}}
        </v-chip> 
        <span>missed cleavages (max)</span>
      </v-tooltip>
      <v-tooltip top>
        <v-chip label outline color="primary" slot="activator">
            <v-icon left>toll</v-icon>{{maxOxidizedMethionine}}
        </v-chip>
        <span>oxidized methionine (max)</span>
      </v-tooltip>
    </div>
    Output format
    <v-radio-group v-model="outputFormat" @click.native="changeTaskFromat" column>
      <v-radio color='primary'
        v-for="f in outputFormats"
        :key="f"
        :label="getFormatText(f)"
        :value="f"
      ></v-radio>
    </v-radio-group>
  <v-card-actions>
  <v-spacer></v-spacer>
  <v-btn small color="accent" @click="step = 4"><v-icon>chevron_left</v-icon>back</v-btn>
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
import CSV from '@/components/CSV.vue'


export default {
  components: {Upload, Submit, Info, CSV},
  props: ['taskid'],
  data () {
    return {
      step: 1,
      collisionEnergy: 30,
      fromFasta: false,
      protease: "Trypsin",
      maxMissedClevages: 0,
      maxOxidizedMethionine: 1,
      proteases: ["Trypsin", "Chymotrypsin", "Glu-C", "Lys-C"],
      missedClevages: [0, 1, 2],
      oxidizedMethionine: [0, 1, 2],
      outputFormat: 'msp',
      outputFormats: ['msp', 'spectronaut'],
      modelIntensityName: false,
      modelIntensityList: false,
      modelIRTName: false,
      modelIRTList: false,
      tag: false,
    }
  },
  methods: {
    getFormatText(format){
      if(format == 'msp')
        return "NIST .MSP Text Format of individual spectra (Skyline and MSPepSearch compatible)";
      if(format == 'spectronaut')
        return "Generic text (Spectronaut compatible). All fragments are reported.";
    },
    changeTaskFromat: function(){
      this.$store.commit('changeFormat', this.outputFormat);
    },
    toStepTwo: function (){
      this.setOptions();
      this.step = 2;
    },
    toStepThree: function (){
      if(this.fromFasta)
        this.$store.commit('reportUpload', {fileType: 'peptides.csv', isSuccessful: false});
      else
        this.$store.commit('reportUpload', {fileType: 'FASTA', isSuccessful: false});
      this.step = 3;
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
    setOptions: function () {
      if(this.fromFasta) {
        this.$store.commit('changeOptions', { 
          'collisionEnergy': this.collisionEnergy,
          'protease': this.protease,
          'missedCleavages': this.maxMissedClevages,
          'oxidizedMethionine': this.maxOxidizedMethionine 
        });
      }
    },
    setModels: function(response){
      //let modelIRTList;
      //this.modelIntensityList = response.data.models.intensities;
      let modelIntensityLists = response.data.models.intensities;
      let modelIRTLists = response.data.models.iRT;
      this.modelIRTList = [];
      //modelIRTList = response.data.models.iRT;
      this.modelIntensityList = [];
      for (let i = 0; i < modelIntensityLists.length; i++){
        //if(modelIntensityLists[i].enabled.spectral_library)
        this.modelIntensityList.push(modelIntensityLists[i]);
      }
      for (let i = 0; i < modelIRTLists.length; i++){
        //if(modelIRTLists[i].enabled.spectral_library)
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
  computed: {
    doneUploads: function () {
      let uploads = this.$store.state.task.uploads;
      if(this.fromFasta)
        return uploads.FASTA
      else
        return uploads['peptides.csv'];
    }
  },
  mounted() {
    this.$store.commit('changeFormat', this.outputFormat);
    let modelsUrl = '/prosit/api/models.xsjs';
    axios.get(modelsUrl).then(this.setModels);
  }

}
</script>