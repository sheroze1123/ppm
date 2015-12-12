#ifndef SERVER_H
#define SERVER_H

#include <string>

class Server {
  public:
    Server(const short port, const double L, const int N, const int N_p,
               const double *mass);

    Server(const Server&) = delete;

    ~Server();

    void marshal(const double *position);

  private:
    int N_p_;
    int server_sock_;
    int conn_sock_;

    void sendln(const std::string& s);
};

#endif // SERVER_H
