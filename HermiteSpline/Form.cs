using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace HermiteSpline
{
    public partial class Form1 : System.Windows.Forms.Form
    {
        double a = 1.0, b = 7.0, nu=1e-3;
        int linknum = 1, funcnum = 4;
        Spline s;

        Thread evaluating;
        
        Dictionary<TrackBar, Chart> tbtochart = new Dictionary<TrackBar, Chart>();

        public Form1()
        {
            InitializeComponent();
            //this.Width = Screen.PrimaryScreen.Bounds.Size.Width;
            //this.Height = Screen.PrimaryScreen.Bounds.Size.Height;
            //
            //toolTipPos.SetToolTip(this.chart1, "shit1");
            //toolTipPos.SetToolTip(this.chart2, "shit2");
            //toolTipPos.SetToolTip(this.chart3, "shit3");
            //toolTipPos.SetToolTip(this.chart4, "shit4");

            tbtochart.Add(trackBar1, chart1);
            tbtochart.Add(trackBar2, chart2);
            tbtochart.Add(trackBar3, chart3);
            tbtochart.Add(trackBar4, chart4);
        }

        private void textBoxBorderA_Leave(object sender, EventArgs e)
        {
            if (Double.TryParse(((TextBox)sender).Text, out double tmp))
                a = tmp;
            else
                ((TextBox)sender).Text = a.ToString();
        }

        private void textBoxBorderB_Leave(object sender, EventArgs e)
        {
            if (Double.TryParse(((TextBox)sender).Text, out double tmp))
                b = tmp;
            else
                ((TextBox)sender).Text = b.ToString();
        }

        private void comboBoxFunctions_SelectedIndexChanged(object sender, EventArgs e)
        {
            funcnum = ((ComboBox)sender).SelectedIndex;
        }

        private void trackBar_Scroll(object sender, EventArgs e)
        {
            TrackBar tb = (TrackBar)sender;
            Chart chart;
            tbtochart.TryGetValue(tb, out chart);
            chart.ChartAreas[0].AxisX.ScaleView.ZoomReset();
            chart.ChartAreas[0].AxisY.ScaleView.ZoomReset();
            if (tb.Value == 0)
                return;


            double xfrom = chart.ChartAreas[0].AxisX.ScaleView.ViewMinimum,
                xto = chart.ChartAreas[0].AxisX.ScaleView.ViewMaximum,
                yfrom = chart.ChartAreas[0].AxisY.ScaleView.ViewMinimum,
                yto = chart.ChartAreas[0].AxisY.ScaleView.ViewMaximum;
            double xk = (xto - xfrom) * ((double)tb.TickFrequency / (double)(tb.Maximum * 2+1));
            double yk = (yto - yfrom) * ((double)tb.TickFrequency / (double)(tb.Maximum * 2+1));

            chart.ChartAreas[0].AxisX.ScaleView.Zoom(xfrom + tb.Value * xk, xto - tb.Value * xk);
            chart.ChartAreas[0].AxisY.ScaleView.Zoom(yfrom + tb.Value * yk, yto - tb.Value * yk);

            //textBox1.Text = chart.ChartAreas[0].AxisX.ScaleView.Position.ToString();
            //textBox2.Text = chart.ChartAreas[0].AxisX.ScaleView.ViewMaximum.ToString();
            //chart.ChartAreas[0].AxisX.ScaleView.ZoomReset();
            //chart.ChartAreas[0].AxisY.ScaleView.ZoomReset();
            //chart.ChartAreas[0].AxisX.ScaleView.Zoom();


            /*
            PropertyValueChangedEventArgs earg = (PropertyValueChangedEventArgs)e;
            textBox1.Text = earg.OldValue.ToString();
            textBox2.Text = sender.ToString();
            */


            /*
            textBox1.Text = chart1.ChartAreas[0].AxisX.ScaleView.Size.ToString();
            TrackBar tb = ((TrackBar)sender);
            switch (tb.Name)
            {
                case "trackBar1":
                    if (tb.Value != 0) {
                        //chart1.ChartAreas[0].AxisX.ScaleView.;
;
                        
                        //chart1.ChartAreas[0].AxisX.ScaleView.Size = (10 - tb.Value + 1) * 0.1;
                        //chart1.ChartAreas[0].AxisY.ScaleView.Size = (10 - tb.Value + 1) * 0.1; 
                        
                    }
                    else
                    {
                        chart1.ChartAreas[0].AxisX.ScaleView.ZoomReset();
                        chart1.ChartAreas[0].AxisY.ScaleView.ZoomReset();
                    }
                    break;
            }
            */
        }

        private void textBoxX_TextChanged(object sender, EventArgs e)
        {
            if(!Double.TryParse(((TextBox)sender).Text, out double x))
                return;
            if (s == null)
                return;





            richTextBoxOutputs.Text = "";
            richTextBoxOutputs.Text += "f(x)    == " + s.OriginFunc(x).ToString("G13") + "\n";
            richTextBoxOutputs.Text += "S(A,x)  == " + s.Eval(x).ToString("G13") + "\n";
            richTextBoxOutputs.Text += "f`(x)   == " + s.OriginDer(x).ToString("G13") + "\n";
            richTextBoxOutputs.Text += "S`(A,x) == " + s.EvalDer(x).ToString("G13") + "\n";

        }

        private void EvalAll()
        {
            /*
            double num = Double.Parse(textBox1.Text);
            double[] nums = new double[42];
            for (int i = 0; i < 42; i++)
                nums[i] = num;
            textBox2.Text = shit(nums, 42).ToString();
            */
            /*
            Int64 ptr = alloc();
            textBox1.Text = ptr.ToString();
            textBox2.Text = fre(ptr).ToString();
            */
            
        }

        private void chart_MouseMove(object sender, MouseEventArgs e)
        {
            Chart chart = (Chart)sender;
            try
            {
                double x = chart.ChartAreas[0].AxisX.PixelPositionToValue((double)chart.PointToClient(MousePosition).X);
                double y = chart.ChartAreas[0].AxisY.PixelPositionToValue((double)chart.PointToClient(MousePosition).Y);
                ;
                toolTipPos.SetToolTip(chart, '{' + x.ToString("G5") + ", " + y.ToString("G5") + '}');
            }
            catch
            {
                return;
            }
        }

        private void textBoxNu_Leave(object sender, EventArgs e)
        {
            nu = Double.Parse(((TextBox)sender).Text);
            if (nu == 0.0)
                nu = Double.PositiveInfinity;
        }

        private void buttonEvalAll_Click(object sender, EventArgs e)
        {

            this.Cursor = Cursors.WaitCursor;
            //evaluating = new Thread(EvalAll);
            //evaluating.Start();
            //evaluating.Join();
            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            try
            {
                s = new Spline(funcnum, linknum, a, b, nu);
            }
            catch (Exception exc)
            {
                MessageBox.Show(exc.Message, "Помилка", MessageBoxButtons.OK, MessageBoxIcon.Error);
                this.Cursor = Cursors.Default;
                return;
            }
            sw.Stop();
            richTextBoxOutputs.Text = sw.Elapsed.ToString() + '\n';
            this.Cursor = Cursors.AppStarting;

            dataGridViewParams.Columns["Column4"].Visible = s.param_count == 5;
            dataGridViewParams.Rows.Clear();
            dataGridViewParams.Rows.Add(s.link_count);
            for (int i = 0; i < s.link_count; i++)
            {
                dataGridViewParams.Rows[i].Cells["Num"].Value = (i+1).ToString();
                dataGridViewParams.Rows[i].Cells["Borders"].Value = "[" + s.X[i].ToString("0.00") + ", " + s.X[i + 1].ToString("0.00") + "]";
                dataGridViewParams.Rows[i].Cells["MaxError"].Value = s.MaxError(s.X[i], s.X[i + 1]).ToString("0.00000000E+00");
                //MaxError
                for (int j = 0; j < s.param_count; j++)
                {
                    dataGridViewParams.Rows[i].Cells[j + 3].Value = s.A[i, j].ToString("0.0000E+00");
                }
            }

            chart1.Series[0].Points.Clear();
            chart1.Series[1].Points.Clear();
            chart2.Series[0].Points.Clear();
            chart2.Series[1].Points.Clear();
            chart3.Series[0].Points.Clear();
            chart3.Series[1].Points.Clear();
            chart4.Series[0].Points.Clear();
            chart4.Series[1].Points.Clear();


            //Color[] cs = { Color.Red, Color.Green, Color.Blue };

            for (int k = 0; k < s.link_count; k++) { 
                double step = (s.X[k+1] - s.X[k]) / (2000);
                for (int i = 0; i <= 2000; i++)
                {
                    double x = s.X[k] + i * step;
                    double y1 = s.OriginFunc(x);
                    double y2 = s.Eval(x);
                    double dy1 = s.OriginDer(x);
                    double dy2 = s.EvalDer(x);
                    chart1.Series[0].Points.AddXY(x, y1);
                    chart1.Series[1].Points.AddXY(x, y2);
                    chart2.Series[0].Points.AddXY(x, y1 - y2);
                    chart3.Series[0].Points.AddXY(x, dy1);
                    chart3.Series[1].Points.AddXY(x, dy2);
                    chart4.Series[0].Points.AddXY(x, dy1 - dy2);
                }
            }
            
            foreach (double x in s.X)
            {
                double y1 = chart2.Series[0].Points.FindByValue(x, "X").YValues[0];
                chart2.Series[1].Points.AddXY(x,y1);
                double y2 = chart4.Series[0].Points.FindByValue(x, "X").YValues[0];
                chart4.Series[1].Points.AddXY(x, y2);
            }

            this.Cursor = Cursors.Default;

            //textBoxX_TextChanged(textBoxX, null);
            //textBox1.Text = chart1.ChartAreas[0].AxisX.ScaleView.Size.ToString();
            //textBox2.Text = chart1.ChartAreas[0].AxisY.ScaleView.Size.ToString();

            //buttonEvalAll.Text = "Побудувати сплайн";
            //buttonEvalAll.BackColor = SystemColors.Control;












            /*
            if (((Button)sender).Text != "Break")
            {
                evaluating = new Thread(EvalAll);
                evaluating.Start();
                ((Button)sender).Text = "Break";
                ((Button)sender).BackColor = Color.Red;
            }
            else
            {
                evaluating.Abort();
                buttonEvalAll.Text = "Побудувати сплайн";
                buttonEvalAll.BackColor = SystemColors.Control;
            }
            */
        }

        /*
        private void numericUpDownNu_ValueChanged(object sender, EventArgs e)
        {
            int newval = (int)(((NumericUpDown)sender).Value);

            //nu = Math.Pow(nu, updown_val - newval);
            if (newval > updown_val)
                nu /= 10;
            else
                nu *= 10;


            textBoxNu.Text = nu.ToString();
            updown_val = newval;
        }
        */

        private void pictureBox_Click(object sender, EventArgs e)
        {
            switch (((PictureBox)sender).Name)
            {
                case "pictureBoxSE4":
                    radioButtonSE4.Checked = true;
                    linknum = 1;
                    break;
                case "pictureBoxSE5":
                    radioButtonSE5.Checked = true;
                    linknum = 2;
                    break;
                case "pictureBoxPN4":
                    radioButtonPN4.Checked = true;
                    linknum = 3;
                    break;
                case "pictureBoxPN5":
                    radioButtonPN5.Checked = true;
                    linknum = 4;
                    break;
            }
        }



    }
}
