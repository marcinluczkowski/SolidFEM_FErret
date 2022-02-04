using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SolidFEM.Classes
{
    class testClass
    {
        public string Name { get; private set; }
        public int ID { get; set; }

        public testClass()
        {
            // empty


        }

        public testClass(string name, int id)
        {
            Name = name;
            ID = id;
        }

        public testClass(testClass copyClass)
        {
            Name = copyClass.Name;
            ID = copyClass.ID;
        }
        

        private void SetName(string name)
        {
            Name = name;
        }
    }
}
