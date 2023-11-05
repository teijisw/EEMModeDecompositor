//EEG MODE DECOMPOSITOR v2.7N
// By Teiji Sawa, MD, PhD
// For eeg_bis.TSV data file
//Nov 4, 2023;
String subscript_1 = "v2.7N by Teiji SAWA, Nov-3-23.";
String subscript_2 = "Anesthesiology, Kyoto Prefectural University of Medicine.";
// ---------------------------------------------
//Swing subwindows
// ---------------------------------------------
int saveFrame = 1;
// ---------------------------------------------
// WMD //--Meyer filter bank
//      boundaries[0] = 4.0D;
//      boundaries[1] = 8.0D;
//      boundaries[2] = 14.0D;
//      boundaries[3] = 20.0D;
//      boundaries[4] = 30.0D;
//---END

double Frequency = 128;

//------------------------------------------------
//int Comp_30 = 0  or Comp_30 = 1 > complementation30 min ;
int Comp_30 = 0;
//------------------------------------------------

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.text.BadLocationException;

StandardOutTest standardOut;
JSelect selector;

int color_trend_size = 256;
int Mode_Flag = 2;
int Start_Flag = 0;
int Stop_Flag = 0;
int Final_Select_Flag = 0;
int loop_cnt1 = 0;
String D_Mode;
String Selected_IMFK;
String Selected_IMF1;
String Selected_IMF2;
int choice=0;
int loop_cnt_for_saveFrame = 0;

String input_datafile;
//INPUT DATA FILE NAME------------------
//String input_datafile = "eeg_bis_20201009_sev_EME30min.tsv";   //Sev-1
//String input_datafile = "eeg_bis_20201016_sev_EVE30min.tsv";   //Sev-2
//String input_datafile = "eeg_bis_20220621_sev_EME30min.tsv";   //Sev-3
//String input_datafile = "eeg_bis_20220627_sev_EME30min.tsv";   //Sev-4
//String input_datafile = "eeg_bis_20220628_1_sev_EME30min.tsv"; //Sev-5
//String input_datafile = "eeg_bis_20220628_2_sev_EME30min.tsv"; //Sev-6
//String input_datafile = "eeg_bis_20220629_sev_EME30min.tsv"; //Sev-7
//String input_datafile = "Sev_8_eeg_bis_30min.tsv"; //Sev-8
//String input_datafile = "Sev_9_eeg_bis_30min.tsv"; //Sev-9
//String input_datafile = "Sev_10_eeg_bis_30min.tsv"; //Sev-10
//String input_datafile = "Sev_11_eeg_bis_30min.tsv"; //Sev-11

double target_freq = 12;
int select_imf;
double select_fq_mean;
double select_fq_stdev;
//--------------------------------------
//Right Graph: Selection for IMF
// Ex)  Gf_R1 = 1, 2, 3, 4, 5, 6, #7(3-4), #8(5-6), 9(all, 1-6), [10(Gf_R1), 11(FFT)], None
int Gf_R1;
//int Gf_R1 = 2;
int Gf_R2;
//int Gf_R1 = 3;
int dual = 0; //spectrogram: 0:single draw, 1: double draw
//--------------------------------------
//Right Graph: Y-axis MAX
// Ex)  Gf_R1_max = 1000
int Gf_E_max = 1000;
int Gf_R1_max = 1000;
//int Gf_E_max;
//int Gf_R1_max;

float limitter = 30;

Emd2 emd_1;
Vmd vmd_1;
Ewt ewt_1;
Ewt2 ewt_2;

//Initial Variables---------------------
int N = 2000;
Complex alpha = new Complex(2000.0, 0);
Complex tau = ZERO;
int K = 6;
int DC = 0;
int init = 1;
double tol = 1e-7;

int T1 = 1024;
//INPUT DATA FILE NAME------------------
public static final Complex ZERO = new Complex(0, 0);
public static final Complex ONE = new Complex(1, 0);
public static final Complex TWO = new Complex(2, 0);
public static final Complex HALF = new Complex(0.5, 0);

FastFourierTransformer fft;
FastFourierTransformer fft1;
FastFourierTransformer fft2;
FastFourierTransformer fft3;
TransformType forward;
TransformType inverse;

PrintWriter output;
PrintWriter output2;
double signal2[];
float t[];
double[][] u_double;


//-------------------------------------------------------------------------------
// For Hilbert calculation and Convolution, "JDSP: a library of digital signal processing tools written in Java"
// https://jdsp.dev
// jdsp-0.5.0 was utilized as a jar file.
// jdsp-0.5.0.jar in "code" folder
// https://github.com/psambit9791/jdsp/tree/master
// @software{sambit_paul_2023_7675362,
//  author       = {Sambit Paul and
//                  Sibo Van Gool},
//  title        = {psambit9791/jdsp: v2.0.1 (February 24, 2023)},
//  month        = feb,
//  year         = 2023,
//  publisher    = {Zenodo},
//  version      = {v2.0.1},
//  doi          = {10.5281/zenodo.7675362},
//  url          = {https://doi.org/10.5281/zenodo.7675362}
// }
//-------------------------------------------------------------------------------
import com.github.psambit9791.jdsp.transform.Hilbert;
//-------------------------------------------------------------------------------
// For FFT and Complex calculation, "Commons Math: The Apache Commons Mathematics Library"
// https://commons.apache.org/proper/commons-math/
// (commons-math3-3.6.1) was utilized as a jar file.
// commons-math3-3.6.1.jar in "code" folder
//-------------------------------------------------------------------------------
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.stat.Frequency;
import org.apache.commons.math3.transform.DftNormalization;

// --------------------------------------------------------------
// EMD: Empirical mode decomposition
// --------------------------------------------------------------
// The implementation for Processing using Java package tryout.emd which was build from the following Java code
// https://stackoverflow.com/questions/10230106/is-there-some-empirical-mode-decomposition-library-in-java
// The above code is from C implementation ( https://code.google.com/p/realtime-emd/ ) and translated it to Java.
// --------------------------------------------------------------
import tryout.emd.Emd;
import tryout.emd.EmdDataImpl;

PImage img, img0, img1, img2, img3, img4, img5, img6, img_anesth_kpum;
PImage imgL0, imgL1, imgL2, imgL3, imgL4, imgL5, imgL6;

FloatTable eeg_data;
float his, his_low;
Graph graph_A, graph_B1, graph_B2, graph_B3, graph_B4, graph_B5, graph_B6, graph_B7, graph_B8,
  graph_C, graph_D, graph_E, graph_F, graph_G, graph_H1, graph_H2,
  graph_H2_1, graph_H2_2, graph_H2_3, graph_H2_4, graph_H2_5, graph_H2_6, graph_H2_7,
  graph_J;

Hilbert hilb[];
EmdDataImpl emd;
double [][] ewt;

int year, month, day, hour, min, sec;
String date_now, time_now, date, time, start_time;
String TitleData_fft, TitleData_hht, TitleData_hht_spec;
int rowCount;
int columnCount;
int current_row=0;

int Hz_64[];
int Fs;
float data_n[];
double x1[];
double x1_r1024[];
float x1_float[];
//double Hamming[];
double Hanning[];
double Blacman[];
//double x1_Hamming[];
double x1_Hanning[];
double x1_Blacman[];

float x1_Hz[];
Complex y[][];
double y_abs[][];
double y_abs_Hz[][];
double y_dB_Hz[][];
float y_abs_Hz_float[][];
float y_dB_Hz_float[][];

double[][][] analytical_signal;
double[][] envelope;
float[][] envelope_float;
double[][] phase;
double[][] frequency;
float[][] frequency_float;
double[][] signal;
float[][] signal_float;
float[][] amp_float;

float [] Freq;
float [] Amp;
float Freq_adj[][];
float Freq_adj2[][][];
float Pw_adj[][];
float Pw_adj2[][][];
float Pw_adj_max[];
float Pw_adj_min[];

float tp_sum;
float tp_sum_low;
float tp[];
float tp_low[];
float fq_mean_avg;
float fq_mean_low_avg;
float fq_mean[];
float fq_mean_low[];
float fq_stdev_avg;
float fq_stdev_low_avg;
float fq_stdev[];
float fq_stdev_low[];
double[] data;
double[][] imfs;
double[][] imfs_temp;
float[][] imfs_float;

PrintWriter output_fft_abs;
PrintWriter output_fft_dB;
PrintWriter output_hht;
PrintWriter output_env;
PrintWriter output_phase;
PrintWriter output_freq;
PrintWriter output_ana_sig;
PrintWriter output_hht_spec_amp;
PrintWriter output_hht_spec_freq;

PFont plotFont;

///////////////////////////////////
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import java.io.InputStreamReader;
import java.io.FileReader;


void setup() {
  //size(1200, 800, P3D);
  size(1600, 900);

  loadImages();

  surface.setResizable(true);
  surface.setSize(1600, 900);
  //noLoop();

//------------------------------------------------
if (Comp_30 == 1) {
  color_trend_size = 232;
}
else {
  color_trend_size = 256;
}
//------------------------------------------------

  smooth();
  plotFont = createFont("SansSerif", 20);
  textFont(plotFont);

  frameRate(1);

  //BEGIN: Java Stanard Out
  standardOut = new StandardOutTest();
  standardOut.setTitle("Console");
  //END: Java Stanard Out

  selector = new JSelect();
  selector.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
  selector.setBounds(10, 10, 540, 210);
  selector.setTitle("Decomposition Mode");
  selector.setVisible(true);

  frameRate(1);
  smooth();

  //Write out the data
  year = year();
  month = month();
  day = day();
  hour = hour();
  min = minute();
  sec = second();

  date_now = year + "_" + nf(month, 2) + "_" + nf(day, 2) + "_" + hour + "-" + nf(min, 2) + "-" + nf(sec, 2);
  time_now = hour + "-" + nf(min, 2) + "-" + nf(sec, 2);
  start_time = hour + ":" + nf(min, 2) + ":" + nf(sec, 2);

  //Console 1
  println("Date&Time:", date_now);
  
  //data_writer();
  //emd_1 = new Emd2(x1, K);
  //vmd_1 = new Vmd(x1, alpha, tau, K, DC, init, tol, N);
  //ewt_1 = new Ewt(x1, K);
  //ewt_2 = new Ewt2(x1, K);
}

void preset() {

  //println("PRESET");
}

void default_set() {
  
  emd_1 = new Emd2(x1, K);
  vmd_1 = new Vmd(x1, alpha, tau, K, DC, init, tol, N);
  ewt_1 = new Ewt(x1, K);
  ewt_2 = new Ewt2(x1, K);

  println("READ DEFAULT_SET");
  Freq_adj = new float [11][2048];
  Freq_adj2 = new float [11][color_trend_size][2048];
  Pw_adj = new float [11][2048];
  Pw_adj2 = new float [11][color_trend_size][2048];

  array_init_float(Freq_adj[Gf_R1-1]);
  array_init_2d_float(Freq_adj2[Gf_R1-1]);
  array_init_float(Pw_adj[Gf_R1-1]);
  array_init_2d_float(Pw_adj2[Gf_R1-1]);

  Pw_adj_max = new float [11];
  Pw_adj_min = new float [11];

  data_n = new float [1024];
  Freq = new float [1024];
  Amp = new  float [1024];

  imfs = new double [6][2*T1];
  imfs_temp = new double [6][2*T1];

  x1 = new double [1040];
  x1_float = new float [1024];
  x1_Hz = new float [1040];
  x1_r1024 = new double [1024];

  //Hamming = new double [1024];
  Hanning = new double [1024];
  Blacman = new double [1024];
  //x1_Hamming = new double [1024];
  x1_Hanning = new double [1024];
  x1_Blacman = new double [1024];

  y = new Complex [3][1024];
  y_abs = new double [3][1024];
  y_abs_Hz = new double [3][256];
  y_dB_Hz = new double [3][256];
  y_abs_Hz_float = new float [3][256];
  y_dB_Hz_float = new float [3][256];

  //imfs_float = new float[K][1024];
  hilb = new Hilbert[9];
  analytical_signal = new double[9][1024][2];
  frequency = new double [9][1024];
  frequency_float = new float [9][1024];
  phase = new double [9][1024];
  envelope = new double [9][1024];
  envelope_float = new float [9][1024];
  signal = new double [9][1024];
   
  signal_float = new float [9][1024];
  amp_float = new float [9][1024];
  fq_mean = new float [9];
  fq_mean_low = new float [9];
  fq_stdev = new float [9];
  fq_stdev_low = new float [9];
  tp = new float [9];
  tp_low = new float [9];
  
  array_init_2d_double(signal);
  array_init_2d_float(signal_float); 

  data = new double [1025];
  Hz_64 = new int[64];

  imfs_float = new float[7][1024];
  x1 = new double [1040];

  t = new float [2*T1];

  for (int i = 0; i < T1; i++) {
    t[i] = (1/float(T1) * (i + 1));
  }

  for (int i=0; i<64; i++) {
    Hz_64[i] = i;
  }

  eeg_data = new FloatTable("data/" + input_datafile);
  String subscript_2 = input_datafile;

  rowCount = eeg_data.getRowCount();
  columnCount = eeg_data.getColumnCount();
}

void imfs_temp(double [][] _imfs) {
  for (int k=0; k < K; k++) {
    for (int i=0; i < int(T1); i++) {
      imfs[k][i] = _imfs[k][i];
      //u_n1_double[i] = u_double[k-1][i];
    }
  }
}

void imfs(double [][] _imfs) {
  for (int k=0; k < K; k++) {
    for (int i=0; i < int(T1); i++) {
      imfs_float[k][i] = (float) _imfs[k][i];
      //u_n1_double[i] = u_double[k-1][i];
    }
  }
}


void draw() {
  try {

    background(224);
    colorMode(HSB);
    fill(0);
    rectMode(CORNERS);
    noStroke();

    Stop_Flag = selector.getStop_Flag();
    Start_Flag = selector.getStart_Flag();
    D_Mode = selector.getdMode();
    Selected_IMFK = selector.getImfK();
    K = int(Selected_IMFK);
    Selected_IMF1 = selector.getImf1();
    Selected_IMF2 = selector.getImf2();
    selected_imf1();
    selected_imf2();
    Gf_E_max = int(selector.get_Y_Max1());
    Gf_R1_max = int(selector.get_Y_Max2());

    if (Stop_Flag == 1) {
      exit();
    } else {

      if (choice == 1) {

        if (current_row < rowCount-64) {

          get_eegData();

          println("x1.length=:", x1.length);
          
          println("D_Mode=", D_Mode);

          switch(D_Mode) {
          case "EMD" :
            imfs_temp = emd_1.emd(x1);
            break;
          case "VMD" :
            imfs_temp = vmd_1.vmd(x1);
            break;
          case "EWT" :
            imfs_temp = ewt_1.ewt(x1_r1024);
            break;
          case "WMD" :
            imfs_temp = ewt_2.ewt(x1_r1024);
            break;
          default:
            imfs_temp = emd_1.emd(x1);
            break;
          }
          
          println("imfs_temp[0].length=:", imfs_temp[0].length);
          imfs_temp(imfs_temp);
          println("imfs[0]=:", imfs[0][0]);
          imfs(imfs);
          println("imfs[0]=:", imfs[0][0]);
          
          hht(imfs);

          draw_graphs();
          power_spectrum1(11, 410);
          power_spectrum2(10, 410);

         switch(K) {
           case 1 :
             power_spectrum3(1, 730);
             power_spectrum3(9, 730);
             break;
           case 2 :
             power_spectrum3(1, 730);
             power_spectrum3(2, 730);
             power_spectrum3(9, 730);
             break;
           case 3 :
             power_spectrum3(1, 730);
             power_spectrum3(2, 730);
             power_spectrum3(3, 730);
             power_spectrum3(9, 730);
             break;
           case 4 :
             power_spectrum3(1, 730);
             power_spectrum3(2, 730);
             power_spectrum3(3, 730);
             power_spectrum3(4, 730);
             power_spectrum3(9, 730);
             break;
           case 5 :
             power_spectrum3(1, 730);
             power_spectrum3(2, 730);
             power_spectrum3(3, 730);
             power_spectrum3(4, 730);
             power_spectrum3(5, 730);
             power_spectrum3(9, 730);
             break;
           case 6 :
             power_spectrum3(1, 730);
             power_spectrum3(2, 730);
             power_spectrum3(3, 730);
             power_spectrum3(4, 730);
             power_spectrum3(5, 730);
             power_spectrum3(6, 730);
             power_spectrum3(9, 730);
             break;
           default :
             power_spectrum3(1, 730);
             power_spectrum3(2, 730);
             power_spectrum3(3, 730);
             power_spectrum3(4, 730);
             power_spectrum3(5, 730);
             power_spectrum3(6, 730);
             power_spectrum3(9, 730);
             break;
          }

          trend_graph();
          write_FFT();
          write_HHT();
          labels();

          saveFrame("user_data/" + date_now + "_" + D_Mode + "_IMF-" + Gf_R1 + "/frames/########.png");
          //loop_cnt_for_saveFrame += 1;
          //if (loop_cnt_for_saveFrame == saveFrame) {
          //    saveFrame("user_data/" + date_now + "_IMF-" + Gf_R1  + "/frames/########.png");
          //    loop_cnt_for_saveFrame = 0;
          //}

          current_row = current_row + 64;
          if (current_row + 64 > rowCount)
          {
            //Console
            println("current_row=:", current_row, "/", rowCount);
            println("END of DATA");

            stop();
          } else {
            //Console
            println("current_row=:", current_row, "/", rowCount);
            //loop();
          }
        } else {
          stop();
        }

        //        loop_cnt2 += 1;
      } else if (choice == 0) {
        title_show();
        //noLoop();
        if (Start_Flag == 1) {
          preset();
          default_set();
          data_writer();
          choice = 1;

          //Console 2
          println("Decomposition_Mode=: ", D_Mode);
          //Console 3
          //Console 4
          println("Input Datafile=: ", input_datafile);
          println("SATARTED");
          println("Dual Mode=: ", dual);
        }
      }
    }
  }   //try
  catch(RuntimeException e) {
  }
}

void get_eegData() {

  for (int j = 0; j < 65; j++) {
    for (int i = 0; i < 16; i++) {
      //x1[j*16+i] = (double) eeg_data.getFloat(rowCount-1-j, i+1);
      x1[j*16+i] = (double) eeg_data.getFloat(current_row + j +1, i+1);
      x1_Hz[j*16+i] = 0.125*(j*16+i);
    }
  }

  for (int i=0; i<1024; i++) {
    data_n[i] = (float) i;
    x1_r1024[i] = x1[i];
    x1_float[i] = (float) x1[i];
    Hanning[i] = 0.5-0.5*cos(2*PI*i/1023);
    x1_Hanning[i] = Hanning[i] * x1[i];
    Blacman[i] = 0.42+0.5*cos(2*PI*i/1023)+ 0.08*cos(4*PI*i/1023);
    x1_Blacman[i] = Blacman[i] * x1[i];
  }

  FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
  y[0] = fft.transform(x1_r1024, TransformType.FORWARD);
  y[1] = fft.transform(x1_Hanning, TransformType.FORWARD);
  y[2] = fft.transform(x1_Blacman, TransformType.FORWARD);

  for (int j=0; j<3; j++) {
    for (int i=0; i<1024; i++) {
      y_abs[j][i] = y[j][i].abs();
    }
    for (int i=0; i<255; i++) {
      y_abs_Hz[j][i] = (y_abs[j][4*i] +  y_abs[j][4*i+1]  +  y_abs[j][4*i+2] +  y_abs[j][4*i+3]);
      y_dB_Hz[j][i] = 20*log(((float) y_abs_Hz[j][i]/0.1))/log(10);
      y_abs_Hz_float[j][i] = (float) y_abs_Hz[j][i];
      y_dB_Hz_float[j][i] = (float) y_dB_Hz[j][i];
    }
  }
}

void hht(double [][]_imfs) {
  for (int j =0; j<6; j++) {
    //signal[j] = u_double[j];
    signal[j] = _imfs[j];
  }
  //signal[6] = sumArray(signal[2], signal[3]);
  //signal[7] = sumArray(signal[3], signal[4]);
  //signal[7] = sumArray(signal[4], signal[5]);
  //signal[7] = sumArray(sumArray(sumArray(signal[1], signal[2]), signal[3]), signal[4]);
  
  switch(K) {
  case 1 :
    //signal[1] = signal[0];
    //signal[2] = signal[0];
    //signal[3] = signal[0];
    //signal[4] = signal[0];
    //signal[5] = signal[0];
    signal[6] = signal[0];
    signal[7] = signal[0];
    signal[8] = signal[0];
    break;
  case 2 :
    //signal[2] = signal[0];
    //signal[3] = signal[0];
    //signal[4] = signal[0];
    //signal[5] = signal[0];
    signal[6] = signal[0];
    signal[7] = signal[0];
    signal[8] = sumArray(signal[0], signal[1]);
    break;
  case 3 :
    //signal[3] = signal[0];
    //signal[4] = signal[0];
    //signal[5] = signal[0];
    signal[6] = signal[0];
    signal[7] = signal[0];
    signal[8] = sumArray(sumArray(signal[0], signal[1]), signal[2]);
    break;
  case 4 :
    //signal[4] = signal[0];
    //signal[5] = signal[0];
    signal[6] = signal[0];
    signal[7] = signal[0];
    signal[8] = sumArray(sumArray(sumArray(signal[0], signal[1]), signal[2]), signal[3]);
    break;
  case 5 :
    //signal[5] = signal[0];
    signal[6] = signal[0];
    signal[7] = signal[0];
    signal[8] = sumArray(sumArray(sumArray(sumArray(signal[0], signal[1]), signal[2]), signal[3]), signal[4]);
    break;
  case 6 :
    signal[6] = signal[0];
    signal[7] = signal[0];
    signal[8] = sumArray(sumArray(sumArray(sumArray(sumArray(signal[0], signal[1]), signal[2]), signal[3]), signal[4]), signal[5]);
    break;
  default:
    signal[6] = signal[0];
    signal[7] = signal[0];
    signal[8] = sumArray(sumArray(sumArray(sumArray(sumArray(signal[0], signal[1]), signal[2]), signal[3]), signal[4]), signal[5]);
    break;
  }
  
  for (int j =0; j<9; j++) {
    for (int i=0; i<1024; i++) {
      signal_float[j][i] =  (float) signal[j][i];
    }
  }

  Fs = 128; //Sampling Frequency of the original signal

  for (int j =0; j<9; j++) {
    hilb[j] = new Hilbert(signal[j]);
    hilb[j].hilbertTransform();
    analytical_signal[j] = hilb[j].getOutput();
    envelope[j] = hilb[j].getAmplitudeEnvelope();
    phase[j] = hilb[j].getInstantaneousPhase();
    frequency[j] = hilb[j].getInstantaneousFrequency(Fs);
    for (int i=0; i<1024; i++) {
      envelope_float[j][i] =  (float) envelope[j][i];
      frequency_float[j][i] = (float) frequency[j][i];
    }
    for (int i=0; i<1024; i++) {
      amp_float[j][i] =   envelope_float[j][i] * envelope_float[j][i];
    }

    tp[j] =  arrSqSum(envelope_float[j]);
    tp_low[j] = arrSqSum_low(limitter, frequency_float[j], envelope_float[j]);
    fq_mean[j] = arrAvg(frequency_float[j]);
    fq_mean_low[j] = arrAvg_low(frequency_float[j], limitter);
    fq_stdev_low[j] = arrStd_low(frequency_float[j], limitter, fq_mean_low[j]);
    fq_stdev[j] = arrStdev(frequency_float[j]);
  }

  //select_imf = Gf_R1;

  select_fq_mean = fq_mean[Gf_R1-1];
  select_fq_stdev = fq_stdev[Gf_R1-1];

  if (tp[Gf_R1 -1]/1000 >= 200) {
    //select_imf=1;
    image(img6, 100, 2);
    image(imgL6, 1280, 130);
  } else if (tp[Gf_R1-1]/1000 >= 100) {
    //select_imf=2;
    image(img5, 100, 2);
    image(imgL5, 1280, 130);
  } else if (tp[Gf_R1-1]/1000 >= 50) {
    //select_imf=3;
    image(img4, 100, 2);
    image(imgL4, 1280, 130);
  } else if (tp[Gf_R1-1]/1000 >= 25) {
    //select_imf=4;
    image(img3, 100, 2);
    image(imgL3, 1280, 130);
  } else if (tp[Gf_R1-1]/1000 >= 10) {
    //select_imf=5;
    image(img2, 100, 2);
    image(imgL2, 1280, 130);
  } else {
    //select_imf=6;
    image(img1, 100, 2);
    image(imgL1, 1280, 130);
  }
}

void write_FFT() {
  year = year();
  month = month();
  day = day();
  sec = second();
  min = minute();
  hour = hour();
  date =  year + ":" + nf(month, 2) + ":" + nf(day, 2);
  time = nf(hour, 2) + ":" + nf(min, 2) + ":" + nf(sec, 2);

  output_fft_abs.print("ch1:" + TAB + time);
  for (int i=0; i < 128; i++) {
    output_fft_abs.print(TAB);
    output_fft_abs.print(y_abs_Hz[2][i]);
  }
  output_fft_abs.println("");
  output_fft_abs.flush();

  output_fft_dB.print("ch1:" + TAB + time);
  for (int i=0; i < 128; i++) {
    output_fft_dB.print(TAB);
    output_fft_dB.print(y_dB_Hz[2][i]);
  }
  output_fft_dB.println("");
  output_fft_dB.flush();
}

void write_HHT() {
  year = year();
  month = month();
  day = day();
  sec = second();
  min = minute();
  hour = hour();
  date =  year + ":" + nf(month, 2) + ":" + nf(day, 2);
  time = nf(hour, 2) + ":" + nf(min, 2) + ":" + nf(sec, 2);

  switch(Gf_R1) {
  case 1 :
    Freq = frequency_float[0];
    Amp  = envelope_float[0];
    tp_sum = tp[0];
    tp_sum_low = tp_low[0];
    fq_mean_avg = fq_mean[0];
    fq_stdev_avg = fq_stdev[0];
    fq_mean_low_avg = fq_mean_low[0];
    fq_stdev_low_avg = fq_stdev_low[0];
    break;
  case 2 :
    Freq = frequency_float[1];
    Amp  = envelope_float[1];
    tp_sum = tp[1];
    tp_sum_low = tp_low[1];
    fq_mean_avg = fq_mean[1];
    fq_stdev_avg = fq_stdev[1];
    fq_mean_low_avg = fq_mean_low[1];
    fq_stdev_low_avg = fq_stdev_low[1];
    break;
  case 3 :
    Freq = frequency_float[2];
    Amp  = envelope_float[2];
    tp_sum = tp[2];
    tp_sum_low = tp_low[2];
    fq_mean_avg = fq_mean[2];
    fq_stdev_avg = fq_stdev[2];
    fq_mean_low_avg = fq_mean_low[2];
    fq_stdev_low_avg = fq_stdev_low[2];
    break;
  case 4 :
    Freq = frequency_float[3];
    Amp  = envelope_float[3];
    tp_sum = tp[3];
    tp_sum_low = tp_low[3];
    fq_mean_avg = fq_mean[3];
    fq_stdev_avg = fq_stdev[3];
    fq_mean_low_avg = fq_mean_low[3];
    fq_stdev_low_avg = fq_stdev_low[3];
    break;
  case 5 :
    Freq = frequency_float[4];
    Amp  = envelope_float[4];
    tp_sum = tp[4];
    tp_sum_low = tp_low[4];
    fq_mean_avg = fq_mean[4];
    fq_stdev_avg = fq_stdev[4];
    fq_mean_low_avg = fq_mean_low[4];
    fq_stdev_low_avg = fq_stdev_low[4];
    break;
  case 6 :
    Freq = frequency_float[5];
    Amp  = envelope_float[5];
    tp_sum = tp[5];
    tp_sum_low = tp_low[5];
    fq_mean_avg = fq_mean[5];
    fq_stdev_avg = fq_stdev[5];
    fq_mean_low_avg = fq_mean_low[5];
    fq_stdev_low_avg = fq_stdev_low[5];
    break;

  case 9 :
    Freq = frequency_float[8];
    Amp  = envelope_float[8];
    tp_sum = tp[8];
    tp_sum_low = tp_low[8];
    fq_mean_avg = fq_mean[8];
    fq_stdev_avg = fq_stdev[8];
    fq_mean_low_avg = fq_mean_low[8];
    fq_stdev_low_avg = fq_stdev_low[8];
    break;

  default:
    Freq = frequency_float[0];
    Amp  = envelope_float[0];
    tp_sum = tp[0];
    tp_sum_low = tp_low[0];
    fq_mean_avg = fq_mean[0];
    fq_stdev_avg = fq_stdev[0];
    fq_mean_low_avg = fq_mean_low[0];
    fq_stdev_low_avg = fq_stdev_low[0];

    break;
  }

  his = 2.8061*(float)fq_mean_avg+15.589;
  his_low = 2.8061*(float)fq_mean_low_avg+15.589;

  output_hht_spec_freq.print("ch1:" + TAB + time);
  for (int i=0; i < 1024; i++) {
    output_hht_spec_freq.print(TAB);
    output_hht_spec_freq.print(Freq[i]);
  }
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(fq_mean_avg);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(fq_mean[0]);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(fq_mean[1]);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(fq_mean[2]);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(fq_mean[3]);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(fq_mean[4]);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(fq_mean[5]);

  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(fq_mean_low_avg);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(his);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(his_low);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(fq_stdev_avg);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(fq_stdev_low_avg);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(tp_sum);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(tp[0]);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(tp[1]);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(tp[2]);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(tp[3]);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(tp[4]);
  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(tp[5]);

  output_hht_spec_freq.print(TAB);
  output_hht_spec_freq.print(tp_sum_low);
  output_hht_spec_freq.println();
  output_hht_spec_freq.flush();

  output_hht_spec_amp.print("ch1:" + TAB + time);
  for (int i=0; i < 1024; i++) {
    output_hht_spec_amp.print(TAB);
    output_hht_spec_amp.print(Amp[i]);
  }
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(fq_mean_avg);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(fq_mean[0]);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(fq_mean[1]);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(fq_mean[2]);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(fq_mean[3]);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(fq_mean[4]);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(fq_mean[5]);

  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(fq_mean_low_avg);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(his);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(his_low);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(fq_stdev_avg);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(fq_stdev_low_avg);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(tp_sum);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(tp[0]);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(tp[1]);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(tp[2]);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(tp[3]);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(tp[4]);
  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(tp[5]);

  output_hht_spec_amp.print(TAB);
  output_hht_spec_amp.print(tp_sum_low);
  output_hht_spec_amp.println();
  output_hht_spec_amp.flush();
}

void data_writer() {
  year = year();
  month = month();
  day = day();
  hour = hour();
  min = minute();
  sec = second();

  date_now = year + "_" + nf(month, 2) + "_" + nf(day, 2) + "_" + hour + "-" + nf(min, 2) + "-" + nf(sec, 2);
  time_now = hour + "-" + nf(min, 2) + "-" + nf(sec, 2);
  start_time = hour + ":" + nf(min, 2) + ":" + nf(sec, 2);

  output_fft_abs = createWriter("user_data/" + date_now + "_" + D_Mode + "_IMF-" + Gf_R1 + "/output_fft_abs.tsv");
  output_fft_dB = createWriter("user_data/" + date_now + "_" + D_Mode + "_IMF-" + Gf_R1 + "/output_fft_dB.tsv");
  output_hht_spec_amp= createWriter("user_data/" + date_now + "_" + D_Mode + "_IMF-" + Gf_R1 + "/output_hht_spec_amp.tsv");
  output_hht_spec_freq= createWriter("user_data/" + date_now + "_" + D_Mode + "_IMF-" + Gf_R1 + "/output_hht_spec_freq.tsv");

  TitleData_fft = "Ch" + TAB + "Time";
  for (int i =0; i< 128; i++) {
    TitleData_fft += (TAB +"sp[" + i + "]");
  }

  output_fft_abs.println(TitleData_fft);
  output_fft_abs.flush();

  output_fft_dB.println(TitleData_fft);
  output_fft_dB.flush();

  TitleData_hht_spec = "Ch" + TAB + "Time";
  for (int i =0; i< 1024; i++) {
    TitleData_hht_spec += (TAB +"sp[" + i + "]");
  }
  TitleData_hht_spec += (TAB +"freq_mean_selected");
  TitleData_hht_spec += (TAB +"freq_mean_1");
  TitleData_hht_spec += (TAB +"freq_mean_2");
  TitleData_hht_spec += (TAB +"freq_mean_3");
  TitleData_hht_spec += (TAB +"freq_mean_4");
  TitleData_hht_spec += (TAB +"freq_mean_5");
  TitleData_hht_spec += (TAB +"freq_mean_6");
  TitleData_hht_spec += (TAB +"freq_mean_low");
  TitleData_hht_spec += (TAB +"HIS");
  TitleData_hht_spec += (TAB +"HIS_low");
  TitleData_hht_spec += (TAB +"freq_stdev");
  TitleData_hht_spec += (TAB +"freq_stdev_low");
  TitleData_hht_spec += (TAB +"tp_selected");
  TitleData_hht_spec += (TAB +"tp_1");
  TitleData_hht_spec += (TAB +"tp_2");
  TitleData_hht_spec += (TAB +"tp_3");
  TitleData_hht_spec += (TAB +"tp_4");
  TitleData_hht_spec += (TAB +"tp_5");
  TitleData_hht_spec += (TAB +"tp_6");
  TitleData_hht_spec += (TAB +"tp_low");

  output_hht_spec_freq.println(TitleData_hht_spec);
  output_hht_spec_freq.flush();

  output_hht_spec_amp.println(TitleData_hht_spec);
  output_hht_spec_amp.flush();
}

void power_spectrum1(int Gf, int yaxis) {
  for (int cnt = 0; cnt < color_trend_size; cnt++) {
    for (int i = 0; i < 94; i++) {
      strokeWeight(1);
      noStroke();

      //println("Pw_adj2[Gf-1][cnt][i]=:", Pw_adj2[Gf-1][cnt][i]);

      if (Freq_adj2[Gf-1][cnt][i] <= 0) {
        fill(0, 0, 0);
      } else {
        //fill(floor(128-(int)((Pw_adj2[Gf-1][cnt][i]-Pw_adj_min[Gf-1])/(Pw_adj_max[Gf-1]-Pw_adj_min[Gf-1])*2048)), 255, 255, 255);
        fill(128-(int)Pw_adj2[Gf-1][cnt][i], 255, 255, 200);
        //println("color_value=:", (int)Pw_adj2[Gf-1][cnt][i]);

        //30min adjestment
        if (Comp_30 == 1) {
          rect(Freq_adj2[Gf-1][cnt][i]-1.8, yaxis+cnt*1.1, Freq_adj2[Gf-1][cnt][i]+1.8, yaxis+(cnt*1.1)+1.8);
        }
        else {
          rect(Freq_adj2[Gf-1][cnt][i]-1.8, yaxis+cnt, Freq_adj2[Gf-1][cnt][i]+1.8, yaxis+cnt+1.1);
        }        

      }
    }
  }
}

void power_spectrum2(int Gf, int yaxis) {
  for (int cnt = 0; cnt < color_trend_size; cnt++) {
    for (int i = 0; i < 2048; i++) {
      strokeWeight(1);
      noStroke();

      if (Freq_adj2[Gf-1][cnt][i] <= 0) {
        fill(0, 0, 0);
      } else {
        fill(floor(128-(int)((Pw_adj2[Gf-1][cnt][i]-Pw_adj_min[Gf-1])/(Pw_adj_max[Gf-1]-Pw_adj_min[Gf-1])*2048)), 255, 255, 50);
        
        //30min adjestment
        if (Comp_30 == 1) {
          rect(Freq_adj2[Gf-1][cnt][i], yaxis+cnt*1.1, Freq_adj2[Gf-1][cnt][i]+1, yaxis+cnt*1.1+1);
        }
        else {
          rect(Freq_adj2[Gf-1][cnt][i], yaxis+cnt, Freq_adj2[Gf-1][cnt][i]+1, yaxis+cnt+1);
        }

    }
    }
  }
}

void power_spectrum3(int Gf, int yaxis) {
  for (int cnt = 0; cnt < color_trend_size; cnt++) {
    for (int i = 0; i < 2048; i++) {
      strokeWeight(1);
      noStroke();

      if (Freq_adj2[Gf-1][cnt][i] <= 0) {
        fill(0, 0, 0);
      } else {
        fill(floor(128-(int)((Pw_adj2[Gf-1][cnt][i]-Pw_adj_min[Gf-1])/(Pw_adj_max[Gf-1]-Pw_adj_min[Gf-1])*2048)), 255, 255, 50);

        //30min adjestment
        if (Comp_30 == 1) {
          rect(Freq_adj2[Gf-1][cnt][i], yaxis+cnt*1.1/2, Freq_adj2[Gf-1][cnt][i]+1, yaxis+cnt*1.1/2+0.5);
        }
        else {
          rect(Freq_adj2[Gf-1][cnt][i], yaxis+cnt/2, Freq_adj2[Gf-1][cnt][i]+1, yaxis+cnt/2+0.5);
        }
      }
    }
  }
}

void trend_graph() {
  for (int i = 0; i < 11; i++) {
    array_shift_2d(Freq_adj2[i], Freq_adj[i]);
    array_shift_2d(Pw_adj2[i], Pw_adj[i]);
  }
}

void draw_graphs() {
  //graph_A = new Graph(0, 1023, -100, 100, 128, 100, 64, 50, 100, 50, 1050, 30);
  graph_A = new Graph(0, 1023, -100, 100, 128, 100, 64, 50, 100, 50, 1050, 30);
  graph_A.drawAxisLabels("EEG [ms]", "uV");
  graph_A.xdrawVolumeLabels();
  graph_A.ydrawVolumeLabels();
  graph_A.drawxLabels();
  graph_A.drawyLabels();
  graph_A.drawPoints(data_n, x1_float);

  graph_B1 = new Graph(0, 1023, -100, 100, 128, 100, 64, 50, 100, 90, 1050, 30);
  graph_B1.drawAxisLabels("", "IMF-1");
  graph_B1.xdrawVolumeLabels();
  graph_B1.ydrawVolumeLabels();
  graph_B1.drawxLabels();
  graph_B1.drawyLabels();
  graph_B1.drawPoints(data_n, imfs_float[0]);

  graph_B2 = new Graph(0, 1023, -100, 100, 128, 100, 64, 50, 100, 130, 1050, 30);
  if(K>1){
    graph_B2.drawAxisLabels("", "IMF-2");
  }
  graph_B2.xdrawVolumeLabels();
  graph_B2.ydrawVolumeLabels();
  graph_B2.drawxLabels();
  graph_B2.drawyLabels();
  graph_B2.drawPoints(data_n, imfs_float[1]);

  graph_B3 = new Graph(0, 1023, -100, 100, 128, 100, 64, 50, 100, 170, 1050, 30);
  if(K>2){
    graph_B3.drawAxisLabels("", "IMF-3");
  }
  graph_B3.xdrawVolumeLabels();
  graph_B3.ydrawVolumeLabels();
  graph_B3.drawxLabels();
  graph_B3.drawyLabels();
  graph_B3.drawPoints(data_n, imfs_float[2]);

  graph_B4 = new Graph(0, 1023, -100, 100, 128, 100, 64, 50, 100, 210, 1050, 30);
  if(K>3){
    graph_B4.drawAxisLabels("", "IMF-4");
  }
  graph_B4.xdrawVolumeLabels();
  graph_B4.ydrawVolumeLabels();
  graph_B4.drawxLabels();
  graph_B4.drawyLabels();
  graph_B4.drawPoints(data_n, imfs_float[3]);

  graph_B5 = new Graph(0, 1023, -100, 100, 128, 100, 64, 50, 100, 250, 1050, 30);
  if(K>4){
    graph_B5.drawAxisLabels("", "IMF-5");
  }
  graph_B5.xdrawVolumeLabels();
  graph_B5.ydrawVolumeLabels();
  graph_B5.drawxLabels();
  graph_B5.drawyLabels();
  graph_B5.drawPoints(data_n, imfs_float[4]);

  graph_B6 = new Graph(0, 1023, -100, 100, 128, 100, 64, 50, 100, 290, 1050, 30);
  if(K>5){
    graph_B6.drawAxisLabels("", "IMF-6");
  }
  graph_B6.xdrawVolumeLabels();
  graph_B6.ydrawVolumeLabels();
  graph_B6.drawxLabels();
  graph_B6.drawyLabels();
  graph_B6.drawPoints(data_n, imfs_float[5]);
  
  //graph_B6 = new Graph(0, 1023, -100, 100, 128, 100, 64, 50, 100, 330, 850, 30);
  //graph_B6.drawAxisLabels("", "IMF-3-4");
  //graph_B6.xdrawVolumeLabels();
  //graph_B6.ydrawVolumeLabels();
  //graph_B6.drawxLabels();
  //graph_B6.drawyLabels();
  //graph_B6.drawPoints(data_n, signal_float[6]);

  //graph_B7 = new Graph(0, 1023, -100, 100, 128, 100, 64, 50, 100, 370, 850, 30);
  graph_B7 = new Graph(0, 1023, -100, 100, 128, 100, 64, 50, 100, 330, 1050, 30);
  graph_B7.drawAxisLabels("", "Σ IMF");
  graph_B7.xdrawVolumeLabels();
  graph_B7.ydrawVolumeLabels();
  graph_B7.drawxLabels();
  graph_B7.drawyLabels();
  graph_B7.drawPoints(data_n, signal_float[8]);

  //graph_B8 = new Graph(0, 1023, -100, 100, 128, 100, 64, 50, 100, 410, 850, 30);
  //graph_B8.drawAxisLabels("", "HHT-"+Gf_R1);
  //graph_B8.xdrawVolumeLabels();
  //graph_B8.ydrawVolumeLabels();
  //graph_B8.drawxLabels();
  //graph_B8.drawyLabels();
  //graph_B8.drawPoints(data_n, frequency_float[Gf_R1-1]);

  //graph_C = new Graph(0, 47, 0, 20000, 10, 4000, 5, 2000, 100, 650, 300, 80);
  graph_C = new Graph(0, 47, 0, 20000, 10, 4000, 5, 2000, 100, 570, 300, 100);
  graph_C.drawAxisLabels_c("Power Spectrum [Hz]", "uV^2");
  graph_C.xdrawVolumeLabels();
  graph_C.ydrawVolumeLabels();
  graph_C.drawxLabels();
  graph_C.drawyLabels();
  graph_C.drawPoints_c(x1_Hz, y_abs_Hz_float[0], y_abs_Hz_float[1], y_abs_Hz_float[2]);

  //graph_D = new Graph(0, 47, 0, 192, 10, 50, 5, 10, 480, 650, 300, 100);
  graph_D = new Graph(0, 47, 0, 192, 10, 50, 5, 10, 480, 570, 300, 100);
  graph_D.drawAxisLabels_c("Power Spectrum [Hz]", "dB");
  graph_D.xdrawVolumeLabels();
  graph_D.ydrawVolumeLabels();
  graph_D.drawxLabels();
  graph_D.drawyLabels();
  graph_D.drawPoints_c(x1_Hz, y_dB_Hz_float[0], y_dB_Hz_float[1], y_dB_Hz_float[2]);

  //graph_E = new Graph(0, 47, 0, Gf_E_max, 10, int(Gf_E_max/5), 5, int(Gf_E_max/10), 100, 500, 300, 100);
  graph_E = new Graph(0, 47, 0, Gf_E_max, 10, int(Gf_E_max/5), 5, int(Gf_E_max/10), 100, 420, 300, 100);
  graph_E.drawAxisLabels_c("", "Pw1");
  graph_E.xdrawVolumeLabels();
  graph_E.ydrawVolumeLabels();
  graph_E.drawxLabels();
  graph_E.drawyLabels();
  graph_E.drawPoints_e(frequency_float, amp_float);

  //graph_F = new Graph(0, 47, 0, Gf_R1_max, 10, int(Gf_R1_max/5), 5, int(Gf_R1_max/10), 480, 500, 300, 100);
  graph_F = new Graph(0, 47, 0, Gf_R1_max, 10, int(Gf_R1_max/5), 5, int(Gf_R1_max/10), 480, 420, 300, 100);
  graph_F.drawAxisLabels_c("", "Pw2");
  graph_F.xdrawVolumeLabels();
  graph_F.ydrawVolumeLabels();
  graph_F.drawxLabels();
  graph_F.drawyLabels();
  if (dual == 0) {
    graph_F.drawPoints_f1(frequency_float[Gf_R1-1], amp_float[Gf_R1-1]);
  } else {
    graph_F.drawPoints_f2(frequency_float[Gf_R1-1], amp_float[Gf_R1-1], frequency_float[Gf_R2-1], amp_float[Gf_R2-1]);
  }

  //graph_G = new Graph(0, 47, 0, Gf_E_max, 10, int(Gf_E_max/4), 5, int(Gf_R1_max), 100, 480, 300, 15);
  graph_G = new Graph(0, 47, 0, Gf_E_max, 10, int(Gf_E_max/4), 5, int(Gf_R1_max), 100, 400, 300, 15);
  graph_G.drawxLabels();
  graph_G.drawyLabels();
  graph_G.drawPoints_g(frequency_float[8], amp_float[8]);

  //graph_H1: Line hilbert spectrum
  //graph_H1 = new Graph(0, 64, 0, Gf_R1_max, 10, Gf_R1_max, 5, Gf_R1_max, 480, 480, 300, 15);
  graph_H1 = new Graph(0, 47, 0, Gf_R1_max, 10, Gf_R1_max, 5, Gf_R1_max, 480, 400, 300, 15);
  graph_H1.drawxLabels();
  graph_H1.drawyLabels();
  if (dual ==0) {
    graph_H1.drawPoints_h1(frequency_float[Gf_R1-1], amp_float[Gf_R1-1]);
  } else {
    graph_H1.drawPoints_h2(frequency_float[Gf_R1-1], amp_float[Gf_R1-1], frequency_float[Gf_R2-1], amp_float[Gf_R2-1]);
  }

  //graph_H1: hilbert spectrogram
  //graph_H2 = new Graph(0, 47, 0, Gf_R1_max, 10, Gf_R1_max, 5, Gf_R1_max, 820, 490, 320, 256);
  graph_H2 = new Graph(0, 47, 0, Gf_R1_max, 10, Gf_R1_max, 5, Gf_R1_max, 1220, 410, 320, 256);
  graph_H2.xdrawVolumeLabels();
  graph_H2.drawxLabels();
  graph_H2.drawyLabels();

  if (dual == 1) {
    graph_H2.drawPoints_h4(Gf_R1, Gf_R2, frequency_float[Gf_R1-1], amp_float[Gf_R1-1], frequency_float[Gf_R2-1], amp_float[Gf_R2-1]);
  } else {
    graph_H2.drawPoints_h3(Gf_R1, frequency_float[Gf_R1-1], amp_float[Gf_R1-1]);
  }

  //FFT Spectrogram
  graph_J = new Graph(0, 47, 0, 150, 10, 200, 5, 200, 820, 410, 320, 256);
  graph_J.xdrawVolumeLabels();
  graph_J.drawxLabels();
  graph_J.drawyLabels();
  graph_J.drawPoints_h5(x1_Hz, y_dB_Hz_float[2]);


  //HSPG IMF-1
  graph_H2_1 = new Graph(0, 47, 0, Gf_R1_max, 10, Gf_R1_max, 5, Gf_R1_max, 100, 730, 160, 128);
  graph_H2_1.xdrawVolumeLabels();
  graph_H2_1.drawxLabels();
  graph_H2_1.drawyLabels();
  graph_H2_1.drawPoints_h3_1(1, frequency_float[0], amp_float[0]);

  //HSPG IMF-2
  graph_H2_2 = new Graph(0, 47, 0, Gf_R1_max, 10, Gf_R1_max, 5, Gf_R1_max, 300, 730, 160, 128);
  graph_H2_2.xdrawVolumeLabels();
  graph_H2_2.drawxLabels();
  graph_H2_2.drawyLabels();
  graph_H2_2.drawPoints_h3_1(2, frequency_float[1], amp_float[1]);

  //HSPG IMF-3
  graph_H2_3 = new Graph(0, 47, 0, Gf_R1_max, 10, Gf_R1_max, 5, Gf_R1_max, 500, 730, 160, 128);
  graph_H2_3.xdrawVolumeLabels();
  graph_H2_3.drawxLabels();
  graph_H2_3.drawyLabels();
  graph_H2_3.drawPoints_h3_1(3, frequency_float[2], amp_float[2]);

  //HSPG IMF-4
  graph_H2_4 = new Graph(0, 47, 0, Gf_R1_max, 10, Gf_R1_max, 5, Gf_R1_max, 700, 730, 160, 128);
  graph_H2_4.xdrawVolumeLabels();
  graph_H2_4.drawxLabels();
  graph_H2_4.drawyLabels();
  graph_H2_4.drawPoints_h3_1(4, frequency_float[3], amp_float[3]);

  //HSPG IMF-5
  graph_H2_5 = new Graph(0, 47, 0, Gf_R1_max, 10, Gf_R1_max, 5, Gf_R1_max, 900, 730, 160, 128);
  graph_H2_5.xdrawVolumeLabels();
  graph_H2_5.drawxLabels();
  graph_H2_5.drawyLabels();
  graph_H2_5.drawPoints_h3_1(5, frequency_float[4], amp_float[4]);

  //HSPG IMF-6
  graph_H2_6 = new Graph(0, 47, 0, Gf_R1_max, 10, Gf_R1_max, 5, Gf_R1_max, 1100, 730, 160, 128);
  graph_H2_6.xdrawVolumeLabels();
  graph_H2_6.drawxLabels();
  graph_H2_6.drawyLabels();
  graph_H2_6.drawPoints_h3_1(6, frequency_float[5], amp_float[5]);

  //HSPG IMF-all
  graph_H2_7 = new Graph(0, 47, 0, Gf_R1_max, 10, Gf_R1_max, 5, Gf_R1_max, 1300, 730, 160, 128);
  graph_H2_7.xdrawVolumeLabels();
  graph_H2_7.drawxLabels();
  graph_H2_7.drawyLabels();
  graph_H2_7.drawPoints_h3_1(9, frequency_float[8], amp_float[8]);
}

void array_shift_2d(float[][] _array1, float[] _array2) {
  for (int j = 1; j<_array1.length; j++) {
    for (int i = 0; i<_array1[j].length; i++) {
      _array1[_array1.length-j][i] = _array1[_array1.length-j-1][i];
      _array1[0][i] = _array2[i];
    }
  }
}

float[][] array_add(float[][] _array1, float[][] _array2) {
  float _array3[][];
  _array3 = new float[_array1.length][_array1[1].length];

  for (int j = 1; j<_array1.length; j++) {
    for (int i = 0; i<_array1[j].length; i++) {
      _array3[j][i] = _array1[j][i] + _array2[j][i];
    }
  }
  return (_array3);
}


void array_init_float(float[] _array) {
  for (int i=0; i<_array.length; i++) {
    _array[i] = 0;
  }
}

void array_init_2d_float(float[][] _array) {
  for (int j=0; j<_array.length; j++) {
    for (int i=0; i<_array[j].length; i++) {
      _array[j][i] = 0;
    }
  }
}

void array_init_double(double[] _array) {
  for (int i=0; i<_array.length; i++) {
    _array[i] = 0;
  }
}

void array_init_2d_double(double[][] _array) {
  for (int j=0; j<_array.length; j++) {
    for (int i=0; i<_array[j].length; i++) {
      _array[j][i] = 0;
    }
  }
}

float arrAvg_low(float[] _array, float limitter) {
  float array_sum = 0;
  int counter = 0;
  for (int i=0; i<_array.length; i++)
  {
    if (_array[i]<= limitter) {
      array_sum += _array[i];
      counter += 1;
    }
  }
  return (array_sum /counter);
}

float arrStd_low(float[] _array, float limitter, float fq_mean_avg) {
  float std_sum = 0;
  int counter = 0;
  for (int i=0; i<_array.length; i++)
  {
    if (_array[i]<= limitter) {
      std_sum += ((_array[i] - fq_mean_avg)  * (_array[i] - fq_mean_avg));
      counter += 1;
    }
  }
  return (sqrt(std_sum) /sqrt(counter));
}

static final float arrSqSum(float... arr) {
  float sum = 0.0;
  for (final float f : arr)  sum += f * f;
  return sum;
}

static final float arrSqSum_low(float limitter, float[] _frequency_float, float... arr) {
  float sum = 0.0;
  for (int i=0; i< arr.length; i++) {
    if (_frequency_float[i] <= limitter) {
      sum += arr[i] * arr[i];
    }
  }
  return sum;
}

static final float arrAvg(float... arr) {
  float average = 0.0;
  for (final float f : arr)  average += f;
  average /= (float)(arr.length);
  return average;
}

static final float arrStdev(float... arr) {
  float average = 0.0;
  float stdev = 0.0;
  for (final float f : arr)  average += f;
  average /= (float)(arr.length);

  for (final float f : arr)  stdev += (f - average)*(f - average);
  stdev = sqrt(stdev/(float)(arr.length));
  return stdev;
}

public static float[] sumArray_float(float[] _array1, float[] _array2) {
  float  sum[];
  sum = new float [_array1.length];
  for (int i = 0; i < _array1.length; i++) {
    sum[i] = (float) _array1[i] + (float) _array2[i];
  }
  return sum;
}

public static double[] sumArray(double[] _array1, double[] _array2) {
  double  sum[];
  sum = new double [_array1.length];
  for (int i = 0; i < _array1.length; i++) {
    sum[i] =  _array1[i] +  _array2[i];
  }
  return sum;
}

float log10 (float x) {
  return (log(x) / log(10));
}

void stop() {
  output_fft_abs.close();
  output_fft_dB.close();
  output_hht_spec_freq.close();
  output_hht_spec_amp.close();
  super.stop();
}

void labels() {
  fill(0);
  textLeading(15);
  textSize(18);
  textAlign(LEFT, CENTER);
  switch(D_Mode) {
  case "EMD" :
    text("Empirical Mode Decomposition:", 155, 12);
    break;
  case "VMD" :
    text("Variational Mode Decomposition:", 155, 12);
    break;
  case "EWT" :
    text("Empirical Wavelet Transform:", 155, 12);
    break;
  case "WMD" :
    text("Wavelet Mode Decomposition:", 155, 12);
    break;
  default:
    text("Empirical Mode Decomposition:", 155, 12);
    break;
  }
  
  textAlign(RIGHT, CENTER);
  textSize(15);
  text("current/total=:" + current_row + "/" + rowCount, 1570, 15);

  textSize(15);
  text("%read=: " + nfc((float)current_row / (float)rowCount * 100, 2) + " %", 1570, 35);

  textAlign(LEFT, CENTER);
  fill(0);
  textSize(15);
  text("HHT", 25, 425);

  fill(0);
  textSize(15);
  text("FFT", 25, 569);

  fill(0);
  textSize(15);
  
  switch(D_Mode) {
  case "EMD" :
    text("EMD", 25, 780);
    break;
  case "VMD" :
    text("VMD", 25, 780);
    break;
  case "EWT" :
    text("EWT", 25, 780);
    break;
  case "WMD" :
    text("WMD", 25, 780);
    break;
  default:
    text("EMD", 25, 780);
    break;
  }
  text("IMFs", 25, 800);

  fill(0);
  textSize(10);
  text("1024", 1150, 360 + 15);
  stroke(0);
  line(1149, 360, 1149, 360+4);

  fill(#1e90ff);
  textSize(12);
  text("IMF-1~6", 43, 405);

  switch(K) {
    case 1 :        
      textSize(12);
      fill(#00ff00);
      text("IMF-1", 120, 427);
      break;

    case 2 :        
      textSize(12);
      fill(#00ff00);
      text("IMF-1", 120, 427);
      fill(#ffff00);
      text("IMF-2", 168, 427);
      break;

    case 3 :             
      textSize(12);
      fill(#00ff00);      
      text("IMF-1", 120, 427);
      fill(#ffff00);
      text("IMF-2", 168, 427);
      fill(#ff0000);
      text("IMF-3", 216, 427);
      break;

    case 4 :             
      textSize(12);
      fill(#00ff00);
      text("IMF-1", 120, 427);
      fill(#ffff00);
      text("IMF-2", 168, 427);
      fill(#ff0000);
      text("IMF-3", 216, 427);
      fill(#00ffff);
      text("IMF-4", 264, 427);
      break;

    case 5 :             
      textSize(12);
      fill(#00ff00);
      text("IMF-1", 120, 427);
      fill(#ffff00);
      text("IMF-2", 168, 427);
      fill(#ff0000);
      text("IMF-3", 216, 427);
      fill(#00ffff);
      text("IMF-4", 264, 427);
      fill(#ff00ff);
      text("IMF-5", 312, 427);
      break;

    case 6 :             
      textSize(12);
      fill(#00ff00);
      text("IMF-1", 120, 427);
      fill(#ffff00);
      text("IMF-2", 168, 427);
      fill(#ff0000);
      text("IMF-3", 216, 427);
      fill(#00ffff);
      text("IMF-4", 264, 427);
      fill(#ff00ff);
      text("IMF-5", 312, 427);
      fill(#ffffff);
      text("IMF-6", 360, 427);
      break;
  }

  fill(#1e90ff);
  textSize(12);
  text("IMF-" + Gf_R1, 443, 405);

  textLeading(15);
  textAlign(LEFT, CENTER);

  fill(#00ff00);
  textSize(12);
  text("Raw", 110, 575);

  fill(#ffff00);
  textSize(12);
  text("Hanning", 180, 575);

  fill(#ff0000);
  textSize(12);
  text("Blackman", 270, 575);

  fill(#00ff00);
  textSize(12);
  text("Raw", 500, 575);

  fill(#ffff00);
  textSize(12);
  text("Hanning", 570, 575);

  fill(#ff0000);
  textSize(12);
  text("Blackman", 660, 575);
  color col1;
  color col2;

  switch(Gf_R1) {
  case 1 :
    col1= #00ff00;
    break;
  case 2 :
    col1= #ffff00;
    break;
  case 3 :
    col1= #ff0000;
    break;
  case 4 :
    col1= #00ffff;
    break;
  case 5 :
    col1= #ff00ff;
    break;
  case 6 :
    col1= #949593;
    break;
  case 7 :
    col1= #ff8c00;
    break;
  case 8 :
    col1= #b22222;
    break;
  case 9 :
    col1= #ff1493;
    break;
  default:
    col1= #00ff00;
    break;
  }

  switch(Gf_R2) {
  case 1 :
    col2= #00ff00;
    break;
  case 2 :
    col2= #ffff00;
    break;
  case 3 :
    col2= #ff0000;
    break;
  case 4 :
    col2= #00ffff;
    break;
  case 5 :
    col2= #ff00ff;
    break;
  case 6 :
    col2= #949593;
    break;
  case 7 :
    col2= #ff8c00;
    break;
  case 8 :
    col2= #b22222;
    break;
  case 9 :
    col2= #ff1493;
    break;
  default:
    col2= #00ff00;
    break;
  }

  if (dual == 1) {
    textSize(14);
    fill(col1);
    text("IMF-" + Gf_R1, 675, 433);
    textSize(14);
    fill(col2);
    text("IMF-" + Gf_R2, 735, 543);
  } else {
    textSize(14);
    fill(col1);
    text("IMF-" + Gf_R1, 735, 433);
  }

  if (dual == 1) {
    textSize(14);
    fill(col1);
    text("IMF-" + Gf_R1, 1440, 398);
    textSize(14);
    fill(col2);
    text("IMF-" + Gf_R2, 1500, 398);
  } else {
    textSize(14);
    fill(col1);
    text("IMF-" + Gf_R1, 1500, 398);
  }

  // Head indicator
  fill(0);
  textSize(20);
  text(D_Mode, 1365, 255);

  if (dual == 1) {
    textSize(20);
    fill(col1);
    text("IMF-" + Gf_R1, 1365, 280);
    textSize(20);
    fill(col2);
    text("IMF-" + Gf_R2, 1365, 300);
  } else {
    textSize(20);
    fill(col1);
    text("IMF-" + Gf_R1, 1365, 280);
  }
  //--------

  fill(#ff0000);
  textSize(12);
  text("Blackman Window", 820, 398);

  fill(0);
  textLeading(15);
  textSize(14);
  textAlign(LEFT, CENTER);
  text("Hilbert Spectrogram:", 1270, 398);

  textSize(15);
  text("Total Power-"+ Gf_R1 + " : "  + ((int) tp_sum/1000) + " x 10^3", 1140, 35);

  if (dual == 1) {
    textSize(20);
    text( "Dual IMF Mode in HSG", 917, 35);
  } else {
    text( "Single IMF Mode in HSG", 917, 35);
  }

  textSize(15);
  text( "Central Freq-" + Gf_R1 + " : "  + nfc((float) select_fq_mean, 1) + " ± " + nfc((float) select_fq_stdev, 1), 1140, 15);

  //textSize(15);
  //text( "± " + nfc((float) select_fq_stdev, 1), 1070, 35);

  //textSize(15);
  //text("Target Freq: " + target_freq + " Hz: IMF-"+ select_imf, 1147, 15);

  fill(0);
  stroke(40);
  strokeWeight(1);

  textAlign(RIGHT, CENTER);
  textSize(12);
  text("min", 810, 390);
  text("0", 811, 410);
  text("5", 811, 450);
  text("10", 811, 490);
  text("15", 811, 530);
  text("20", 811, 570);
  text("25", 811, 610);
  text("30", 811, 650);

  stroke(0);
  line(818, 410, 818, 666);

  line(816, 410, 818, 410);
  line(816, 450, 818, 450);
  line(816, 490, 818, 490);
  line(816, 530, 818, 530);
  line(816, 570, 818, 570);
  line(816, 610, 818, 610);
  line(816, 650, 818, 650);

  stroke(#6e6e6e, 150);
  line(818, 450, 1140, 450);
  line(818, 490, 1140, 490);
  line(818, 530, 1140, 530);
  line(818, 570, 1140, 570);
  line(818, 610, 1140, 610);
  line(818, 650, 1140, 650);

  //-------HHT Spectrogram
  textAlign(RIGHT, CENTER);
  textSize(12);
  text("min", 1210, 390);
  text("0", 1211, 410);
  text("5", 1211, 450);
  text("10", 1211, 490);
  text("15", 1211, 530);
  text("20", 1211, 570);
  text("25", 1211, 610);
  text("30", 1211, 650);

  stroke(0);
  line(1218, 410, 1218, 666);

  line(1216, 410, 1218, 410);
  line(1216, 450, 1218, 450);
  line(1216, 490, 1218, 490);
  line(1216, 530, 1218, 530);
  line(1216, 570, 1218, 570);
  line(1216, 610, 1218, 610);
  line(1216, 650, 1218, 650);

  stroke(#6e6e6e, 150);
  line(1218, 450, 1540, 450);
  line(1218, 490, 1540, 490);
  line(1218, 530, 1540, 530);
  line(1218, 570, 1540, 570);
  line(1218, 610, 1540, 610);
  line(1218, 650, 1540, 650);
  //----------------------------
  textSize(15);
  textAlign(LEFT, CENTER);
  text("Power Spectrogram [Hz]", 900, 700);

  textSize(15);
  textAlign(LEFT, CENTER);
  text("Hilbert Spectrogram [Hz]", 1300, 700);


  switch(K) {
    case 1 : 
      textSize(12);
      text("IMF-1", 104, 720);
      text("Σ IMF", 1304, 720);
      break;
    case 2 : 
      textSize(12);
      text("IMF-1", 104, 720);
      text("IMF-2", 304, 720);
      text("Σ IMF", 1304, 720);
      break;
    case 3 : 
      textSize(12);
      text("IMF-1", 104, 720);
      text("IMF-2", 304, 720);
      text("IMF-3", 504, 720);
      text("Σ IMF", 1304, 720);
      break;
    case 4 : 
      textSize(12);
      text("IMF-1", 104, 720);
      text("IMF-2", 304, 720);
      text("IMF-3", 504, 720);
      text("IMF-4", 704, 720);
      text("Σ IMF", 1304, 720);
      break;
    case 5 : 
      textSize(12);
      text("IMF-1", 104, 720);
      text("IMF-2", 304, 720);
      text("IMF-3", 504, 720);
      text("IMF-4", 704, 720);
      text("IMF-5", 904, 720);
      text("Σ IMF", 1304, 720);
      break;
    case 6 : 
      textSize(12);
      text("IMF-1", 104, 720);
      text("IMF-2", 304, 720);
      text("IMF-3", 504, 720);
      text("IMF-4", 704, 720);
      text("IMF-5", 904, 720);
      text("IMF-6", 1104, 720);
      text("Σ IMF", 1304, 720);
      break;
    default : 
      textSize(12);
      text("IMF-1", 104, 720);
      text("IMF-2", 304, 720);
      text("IMF-3", 504, 720);
      text("IMF-4", 704, 720);
      text("IMF-5", 904, 720);
      text("IMF-6", 1104, 720);
      text("Σ IMF", 1304, 720);
      break;
  }

  if (D_Mode == "WMD") {
    text("δ [0.5-4Hz]", 154, 720);
    text("θ [4-8Hz]", 354, 720);
    text("α [8-14Hz]", 554, 720);
    text("lo-β [14-20Hz]", 754, 720);
    text("hi-β [20-30Hz]", 954, 720);
    text("γ [30-64Hz]", 1154, 720);
    text("  [0.5-64Hz]", 1354, 720);
  }

  textAlign(RIGHT, CENTER);
  textSize(10);
  text("min", 90, 720);
  text("0", 95, 732);
  text("10", 95, 772);
  text("20", 95, 812);
  text("30", 95, 852);

  stroke(0);
  line(99, 730, 99, 854);
  line(97, 730, 99, 730);
  line(97, 770, 99, 770);
  line(97, 810, 99, 810);
  line(97, 850, 99, 850);

  stroke(#6e6e6e, 150);
  line(100, 770, 260, 770);
  line(100, 810, 260, 810);
  line(100, 850, 260, 850);

  textAlign(RIGHT, CENTER);
  textSize(10);
  text("min", 290, 720);
  text("0", 295, 732);
  text("10", 295, 772);
  text("20", 295, 812);
  text("30", 295, 852);

  stroke(0);
  line(299, 730, 299, 854);
  line(297, 730, 299, 730);
  line(297, 770, 299, 770);
  line(297, 810, 299, 810);
  line(297, 850, 299, 850);

  stroke(#6e6e6e, 150);
  line(300, 770, 460, 770);
  line(300, 810, 460, 810);
  line(300, 850, 460, 850);

  textAlign(RIGHT, CENTER);
  textSize(10);
  text("min", 490, 720);
  text("0", 495, 732);
  text("10", 495, 772);
  text("20", 495, 812);
  text("30", 495, 852);

  stroke(0);
  line(499, 730, 499, 854);
  line(497, 730, 499, 730);
  line(497, 770, 499, 770);
  line(497, 810, 499, 810);
  line(497, 850, 499, 850);

  stroke(#6e6e6e, 150);
  line(500, 770, 660, 770);
  line(500, 810, 660, 810);
  line(500, 850, 660, 850);

  textAlign(RIGHT, CENTER);
  textSize(10);
  text("min", 690, 720);
  text("0", 695, 732);
  text("10", 695, 772);
  text("20", 695, 812);
  text("30", 695, 852);

  stroke(0);
  line(699, 730, 699, 854);
  line(697, 730, 699, 730);
  line(697, 770, 699, 770);
  line(697, 810, 699, 810);
  line(697, 850, 699, 850);

  stroke(#6e6e6e, 150);
  line(700, 770, 860, 770);
  line(700, 810, 860, 810);
  line(700, 850, 860, 850);

  textAlign(RIGHT, CENTER);
  textSize(10);
  text("min", 890, 720);
  text("0", 895, 732);
  text("10", 895, 772);
  text("20", 895, 812);
  text("30", 895, 852);

  stroke(0);
  line(899, 730, 899, 854);
  line(897, 730, 899, 730);
  line(897, 770, 899, 770);
  line(897, 810, 899, 810);
  line(897, 850, 899, 850);

  stroke(#6e6e6e, 150);
  line(900, 770, 1060, 770);
  line(900, 810, 1060, 810);
  line(900, 850, 1060, 850);

  textAlign(RIGHT, CENTER);
  textSize(10);
  text("min", 1090, 720);
  text("0", 1095, 732);
  text("10", 1095, 772);
  text("20", 1095, 812);
  text("30", 1095, 852);

  stroke(0);
  line(1099, 730, 1099, 854);
  line(1097, 730, 1099, 730);
  line(1097, 770, 1099, 770);
  line(1097, 810, 1099, 810);
  line(1097, 850, 1099, 850);

  stroke(#6e6e6e, 150);
  line(1100, 770, 1260, 770);
  line(1100, 810, 1260, 810);
  line(1100, 850, 1260, 850);

  textAlign(RIGHT, CENTER);
  textSize(10);
  text("min", 1290, 720);
  text("0", 1295, 732);
  text("10", 1295, 772);
  text("20", 1295, 812);
  text("30", 1295, 852);

  stroke(0);
  line(1299, 730, 1299, 854);
  line(1297, 730, 1299, 730);
  line(1297, 770, 1299, 770);
  line(1297, 810, 1299, 810);
  line(1297, 850, 1299, 850);
  
  stroke(#6e6e6e, 150);
  line(1300, 770, 1460, 770);
  line(1300, 810, 1460, 810);
  line(1300, 850, 1460, 850);
  
  textSize(26);
  textAlign(LEFT, CENTER);
  text("N_of_IMFs=:" + K, 1300, 80);
 
  textSize(12);
  textAlign(LEFT, CENTER);
  text("Hilbert Spectrogram [Hz]", 720, 890);
  
  textSize(15);
  textAlign(LEFT, CENTER);
  text(subscript_1, 155, 35);
  text(subscript_2, 440, 35);
  text("Data: "+input_datafile, 490, 15);
}

void loadImages() {
  img0 = loadImage("./image/Brain_Sig_IMF0.png");
  img1 = loadImage("./image/Brain_Sig_IMF1.png");
  img2 = loadImage("./image/Brain_Sig_IMF2.png");
  img3 = loadImage("./image/Brain_Sig_IMF3.png");
  img4 = loadImage("./image/Brain_Sig_IMF4.png");
  img5 = loadImage("./image/Brain_Sig_IMF5.png");
  img6 = loadImage("./image/Brain_Sig_IMF6.png");
  imgL0 = loadImage("./image/Brain_Sig_L_IMF0.png");
  imgL1 = loadImage("./image/Brain_Sig_L_IMF1.png");
  imgL2 = loadImage("./image/Brain_Sig_L_IMF2.png");
  imgL3 = loadImage("./image/Brain_Sig_L_IMF3.png");
  imgL4 = loadImage("./image/Brain_Sig_L_IMF4.png");
  imgL5 = loadImage("./image/Brain_Sig_L_IMF5.png");
  imgL6 = loadImage("./image/Brain_Sig_L_IMF6.png");
  img = loadImage("./image/Brain_Sig_IMF.png");
  img_anesth_kpum = loadImage("./image/Anesth-KPUM-logo.png");
}

void title_show() {

  image(img, 200, 160);
  image(img_anesth_kpum, 390, 650);
  fill(0);
  textLeading(15);
  textSize(78);
  textAlign(LEFT, CENTER);
  text("EEG Mode Decompositor", 330, 200);

  textSize(48);
  textAlign(LEFT, CENTER);
  text("View Four Mode Decompositions", 450, 300);

  textSize(28);
  textAlign(LEFT, LEFT);
  fill(#0000ff);
  text("1-EMD: Empirical Mode Decomposition", 450, 370);
  text("2-VMD: Variational Mode Decomposition", 450, 410);
  text("3-EWT: Empirical Wavelet Transform", 450, 450);
  text("4-WMD: Wavelet Mode Decomposition: Hz-fixed", 450, 490);
  text("The number of IMFs: 2~6", 450, 530);
  text("IMF selection: IMF-1, IMF-2, IMF-3, IMF-4, IMF-5, IMF-6", 450, 570);
  
  fill(#FF0000);
  text("Select Data File and Click START!!!", 450, 610);
  fill(0);
  textSize(20);
  textAlign(LEFT, CENTER);
  text(subscript_1, 500, 680);
  textAlign(LEFT, CENTER);
  text(subscript_2, 500, 710);
}

void mousePressed() {
  noLoop();
}

void mouseReleased() {
  loop();
}

void selected_imf1() {
  switch(Selected_IMF1) {
  case "IMF-1" :
    Gf_R1=1;
    break;
  case "IMF-2" :
    Gf_R1=2;
    break;
  case "IMF-3" :
    Gf_R1=3;
    break;
  case "IMF-4" :
    Gf_R1=4;
    break;
  case "IMF-5" :
    Gf_R1=5;
    break;
  case "IMF-6" :
    Gf_R1=6;
    break;
  case "IMF-7" :
    Gf_R1=7;
    break;
  case "IMF-8" :
    Gf_R1=8;
    break;
  case "IMF-9" :
    Gf_R1=9;
    break;
  default:
    Gf_R1=1;
    break;
  }
}

void selected_imf2() {
  switch(Selected_IMF2) {
  case "IMF-1" :
    Gf_R2=1;
    dual = 1;
    break;
  case "IMF-2" :
    Gf_R2=2;
    dual = 1;
    break;
  case "IMF-3" :
    Gf_R2=3;
    dual = 1;
    break;
  case "IMF-4" :
    Gf_R2=4;
    dual = 1;
    break;
  case "IMF-5" :
    Gf_R2=5;
    dual = 1;
    break;
  case "IMF-6" :
    Gf_R2=6;
    dual = 1;
    break;
  default:
    Gf_R2=1;
    dual = 0;
    break;
  }
}
