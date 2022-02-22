using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SolidFEM.Classes
{
    public class FEMLogger
    {
        // -- properties --
        public List<string> LogList { get; private set; }
        public System.Diagnostics.Stopwatch Stopwatch { get; set; } // not sure if helpful

        // -- constructors --
        public FEMLogger()
        {
            // empty constructor
            this.LogList = new List<string>();
        }


        // -- methods --
        public void AddInfo(string _info)
        {
            string msg = "Info: ";
            msg += (_info + "\n");
            LogList.Add(msg);
        }

        public void AddError(string _error)
        {
            string msg = "Error: ";
            msg += (_error + "\n");
            LogList.Add(msg);
        }

        public void AddWarning(string _warning)
        {
            string msg = "Warning: ";
            msg += (_warning + "\n");
            LogList.Add(msg);
        }

        public void LogToTextFile(string fileName ="SolidFEM_Comp_log.txt")
        {
            // test if a logger folder already exists; if not, create it. 
            
            string dir = Environment.CurrentDirectory;
            var logDir = Path.Combine(dir, "Loggers");
            try
            {
                if (!Directory.Exists(logDir))
                {
                    Directory.CreateDirectory(Path.Combine(dir, "Loggers"));
                }
            }
            catch (Exception e)
            {
                throw new Exception("Unable to create path for log files.");
            }


            try
            {
                string pathName = Path.Combine(logDir , fileName);
                //string testPath1 = Directory.GetCurrentDirectory();
                //string testPath2 = Directory.GetParent(Environment.CurrentDirectory).Parent.FullName;
                
                if (File.Exists(pathName))
                {
                    File.Delete(pathName);
                }
                // create a new file
                using (FileStream fs =  File.Create(pathName))
                {
                    // add the logger to the file
                    Byte[] title = new UTF8Encoding(true).GetBytes("Log File \n");
                    fs.Write(title, 0, title.Length);

                    foreach (string line in LogList)
                    {
                        Byte[] logLine = new UTF8Encoding(true).GetBytes(line);
                        fs.Write(logLine, 0, line.Length);
                    }

                    fs.Close();
                }
            }
            
            catch (Exception e)
            {
                throw new Exception("Unable to write logger to file!");
            }
        }


        public int GetLength()
        {
            return LogList.Count;
        }


    }
}
