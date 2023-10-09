#ifndef ALGORITHMS_H
#define ALGORITHMS_H
#include <errno.h>     // for WEXITSTATUS etc
#include<strings.h>    // for strerror
#include <iostream>
#include <string>

int SystemCall(std::string cmd, std::string& errmsg);

template <typename Arg, typename... Args>
void Print(std::ostream& out, Arg&& arg, Args&&... args){
    out << std::forward<Arg>(arg);
    using expander = int[];
    (void)expander{0, (void(out << ',' << std::forward<Args>(args)), 0)...};
}

bool CheckPath(std::string p, std::string& t);
bool checkfileexists(std::string);

#endif
