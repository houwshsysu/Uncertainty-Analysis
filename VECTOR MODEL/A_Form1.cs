using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace VECTOR_MODEL
{
    public partial class A_Form1 : Form
    {
        public A_Form1()
        {
            InitializeComponent();
        }

        private string GetDataPath()
        {
            string path = null;
            FolderBrowserDialog OpenTin = new FolderBrowserDialog();
            DialogResult result = OpenTin.ShowDialog();
            if (result == DialogResult.OK)
            {
                path = OpenTin.SelectedPath;
            }
            return path;

        }
    }
}
