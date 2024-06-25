import javax.swing.*;
import java.awt.Dimension;
import java.awt.BorderLayout;
import java.awt.event.*;
import processing.serial.*;
import javax.swing.filechooser.*;
import javax.swing.filechooser.FileFilter;

public class JSelect extends javax.swing.JFrame implements ActionListener {
  JFrame frame = new JFrame("Title");
  JComboBox combo1;
  JComboBox combo2;
  JComboBox combo3;
  JComboBox combo4;
  JComboBox combo5;
  JComboBox combo6;
  JButton button_file_select;
  JButton button_start;
  JButton button_stop;

  JLabel label;
  JLabel label2;
  JLabel label3;
  JTextArea textarea2;
  String data_path;
  File file;
  //File file_tsv;
  // File input_datafile;

  String selected_mode;
  String selected_imfk ="6";
  String selected_dmode ="EMD";
  String Y_Max1 = "1000";
  String Y_Max2 = "1000";
  String selected_imf1 = "IMF-1";
  String selected_imf2 = "IMF-1";
  int selected_mode_flag = 0;
  int start_flag = 0;
  int stop_flag = 0;
  int file_select_flag = 0;
  String selectedFileName;
  JLabel l1, l2, l3, l4, l5;

  JSelect() {
    String[] dmodes = {"EMD", "VMD", "EWT", "WMD"};
    String[] imfk = {"6", "5", "4", "3", "2"};
    String[] imfs1 = {"IMF-1", "IMF-2", "IMF-3", "IMF-4", "IMF-5", "IMF-6"};
    String[] imfs2= {"None", "IMF-1", "IMF-2", "IMF-3", "IMF-4", "IMF-5", "IMF-6"};
    String[] y_max1 = {"1000", "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000", "20000"};
    String[] y_max2 = {"1000", "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000", "20000"};

    //String[] modes = {"FileSelect"};
    button_file_select = new JButton("Data File Select");
    button_start = new JButton("Start");
    button_stop = new JButton("Exit");

    combo1 = new JComboBox(dmodes);
    combo1.setPreferredSize(new Dimension(200, 30));

    combo2 = new JComboBox(imfk);
    combo2.setPreferredSize(new Dimension(100, 30));

    combo3 = new JComboBox(imfs1);
    combo3.setPreferredSize(new Dimension(120, 30));

    combo4 = new JComboBox(imfs2);
    combo4.setPreferredSize(new Dimension(120, 30));

    combo5 = new JComboBox(y_max1);
    combo5.setPreferredSize(new Dimension(100, 30));

    combo6 = new JComboBox(y_max2);
    combo6.setPreferredSize(new Dimension(100, 30));

    //textarea2 = new JTextArea();
    //textarea2.setLineWrap(true);

    //JScrollPane scrollpane = new JScrollPane(textarea2);
    //scrollpane.setPreferredSize(new Dimension(400, 400));

    combo1.addActionListener(this);
    combo2.addActionListener(this);
    combo3.addActionListener(this);
    combo4.addActionListener(this);
    combo5.addActionListener(this);
    combo6.addActionListener(this);
    button_file_select.addActionListener(this);
    button_start.addActionListener(this);
    button_stop.addActionListener(this);

    JPanel p = new JPanel();
    l1 = new JLabel("The number of IMF");
    l2 = new JLabel("1st-IMF in HSG");
    l3 = new JLabel("Dual Mode:2nd-IMF in HSG");
    l4 = new JLabel("Left HS Y-Max1");
    l5 = new JLabel("Center HS Y-Max1");

    p.add(new JLabel("Select D-Mode:"));
    p.add(combo1);
    p.add(l1);
    p.add(combo2);
    p.add(l2);
    p.add(combo3);
    p.add(l3);
    p.add(combo4);
    p.add(l4);
    p.add(combo5);
    p.add(l5);
    p.add(combo6);
    p.add(button_file_select);
    p.add(button_start);
    p.add(button_stop);
    //p.add(scrollpane);

    label = new JLabel();
    JPanel labelPanel = new JPanel();
    labelPanel.add(label);

    getContentPane().add(p, BorderLayout.CENTER);
    getContentPane().add(labelPanel, BorderLayout.PAGE_END);
  }

  public String getMode() {
    return selected_mode;
  }

  public String getImfK() {
    return selected_imfk;
  }

  public String getdMode() {
    return selected_dmode;
  }

  public String getImf1() {
    return selected_imf1;
  }

  public String getImf2() {
    return selected_imf2;
  }

  public int getMode_Flag() {
    return selected_mode_flag;
  }

  public String get_Y_Max1() {
    return Y_Max1;
  }

  public String get_Y_Max2() {
    return Y_Max2;
  }

  public int getStart_Flag() {
    return start_flag;
  }

  public int getStop_Flag() {
    return stop_flag;
  }

  public int getFile_Select_Flag() {
    file_select_flag = 1;
    return file_select_flag;
  }

  public File getFile() {
    return file;
  }

  public void actionPerformed(ActionEvent e) {

    if (e.getSource()==button_stop) {
      println("Stop Button Pushed!");
      stop_flag = 1;
      exit();
    } else if (e.getSource()==button_start) {
      if (file_select_flag == 1) {
        println("Start Button Pushed!");
        start_flag = 1;
      } else {
        println("Choose Data File!");
        JLabel label2 = new JLabel("Please Select Data File!");
        label2.setForeground(Color.RED);
        JOptionPane.showMessageDialog(this, label2);
      }
    } else if (e.getSource()== combo1) {
      selected_dmode = (String)combo1.getSelectedItem();
    } else if (e.getSource()==  combo2) {
      selected_imfk = (String)combo2.getSelectedItem();      
    } else if (e.getSource()==  combo3) {
      selected_imf1 = (String)combo3.getSelectedItem();
      if(selected_imf1 == "IMF-6" && int(selected_imfk) < 6){
        JLabel label3 = new JLabel("Selected IMF over n_IMF!");
        label3.setForeground(Color.RED);
        JOptionPane.showMessageDialog(this, label3);
      }
      if(selected_imf1 == "IMF-5" && int(selected_imfk) < 5){
        JLabel label3 = new JLabel("Selected IMF over n_IMF!");
        label3.setForeground(Color.RED);
        JOptionPane.showMessageDialog(this, label3);
      }   
      if(selected_imf1 == "IMF-4" && int(selected_imfk) < 4){
        JLabel label3 = new JLabel("Selected IMF over n_IMF!");
        label3.setForeground(Color.RED);
        JOptionPane.showMessageDialog(this, label3);
      }   
      if(selected_imf1 == "IMF-3" && int(selected_imfk) < 3){
        JLabel label3 = new JLabel("Selected IMF over n_IMF!");
        label3.setForeground(Color.RED);
        JOptionPane.showMessageDialog(this, label3);
      }   
      if(selected_imf1 == "IMF-2" && int(selected_imfk) < 2){
        JLabel label3 = new JLabel("Selected IMF over n_IMF!");
        label3.setForeground(Color.RED);
        JOptionPane.showMessageDialog(this, label3);
      }   
      
    } else if (e.getSource()==  combo4) {
      selected_imf2 = (String)combo4.getSelectedItem();
      if(selected_imf2 == "IMF-6" && int(selected_imfk) < 6){
        JLabel label3 = new JLabel("Selected IMF over n_IMF!");
        label3.setForeground(Color.RED);
        JOptionPane.showMessageDialog(this, label3);
      }
      if(selected_imf2 == "IMF-5" && int(selected_imfk) < 5){
        JLabel label3 = new JLabel("Selected IMF over n_IMF!");
        label3.setForeground(Color.RED);
        JOptionPane.showMessageDialog(this, label3);
      }   
      if(selected_imf2 == "IMF-4" && int(selected_imfk) < 4){
        JLabel label3 = new JLabel("Selected IMF over n_IMF!");
        label3.setForeground(Color.RED);
        JOptionPane.showMessageDialog(this, label3);
      }   
      if(selected_imf2 == "IMF-3" && int(selected_imfk) < 3){
        JLabel label3 = new JLabel("Selected IMF over n_IMF!");
        label3.setForeground(Color.RED);
        JOptionPane.showMessageDialog(this, label3);
      }   
      if(selected_imf2 == "IMF-2" && int(selected_imfk) < 2){
        JLabel label3 = new JLabel("Selected IMF over n_IMF!");
        label3.setForeground(Color.RED);
        JOptionPane.showMessageDialog(this, label3);
      }   
    
    } else if (e.getSource()==  combo5) {
      Y_Max1 = (String)combo5.getSelectedItem();
    } else if (e.getSource()==  combo6) {
      Y_Max2 = (String)combo6.getSelectedItem();
    } else {

      //if ( (String)combo1.getSelectedItem() == "FileSelect") {
      if (e.getSource()== button_file_select) {
        data_path = dataPath("");
        File dir = new File(data_path + "/../data/");
        JFileChooser filechooser = new JFileChooser(dir);
        //int read_buffer_size = 3408;
        //int loop_cnt2 = 0;
        int selected = filechooser.showOpenDialog(this);

        if (selected == JFileChooser.APPROVE_OPTION) {
          FileFilter filter1 = new FileNameExtensionFilter("dat title", "tsv");
          filechooser.addChoosableFileFilter(filter1);
          file = filechooser.getSelectedFile();
          file_select_flag = 1;
          selected_mode_flag = 2;
        }
        selected_dmode = (String)combo1.getSelectedItem();
        selected_imfk = (String)combo2.getSelectedItem();
        selected_imf1 = (String)combo3.getSelectedItem();
        selected_imf2 = (String)combo4.getSelectedItem();
        println("File Selected!");
        //selected_mode_flag = 2;
        selectedFileName = file.getName();
        label.setText("File: " + selectedFileName );
        //Final_Select_Flag = selector.getFile_Select_Flag();
        //file_tsv = selector.getFile();
        input_datafile = selectedFileName;
        //start_flag = 1;
      }
    }
  }
}
