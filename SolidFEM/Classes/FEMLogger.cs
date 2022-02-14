using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SolidFEM.Classes
{
    public class FEMLogger
    {
        // -- properties --
        public List<string> LogList { get; private set; }

        // -- constructors --
        public FEMLogger()
        {
            // empty constructor
        }


        // -- methods --
        public void AddInfo(string _info)
        {
            LogList.Add(_info);
        }



        public int GetLength()
        {
            return LogList.Count;
        }


    }
}
