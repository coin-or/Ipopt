#include "IpIpoptType.hpp"

namespace Ipopt {

  std::list<IpoptTypeInfo*>& IpoptTypeInfosList()
  {
    static std::list<IpoptTypeInfo*> ipopt_type_infos;
    return ipopt_type_infos;
  }


}
