#include "server.h"

#include <string>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>

Server::Server(const short port, const double L, const int N, const int N_p,
               const double *mass)
    : N_p_(N_p),
      server_sock_(0),
      conn_sock_(0)
{
    // TODO: error checking
    std::string port_string = std::to_string(port);
    struct addrinfo hints;
    struct addrinfo *res;
    struct sockaddr_storage their_addr;
    socklen_t addr_size = sizeof(their_addr);

    memset(&hints, 0, sizeof(hints)); // make sure the struct is empty
    hints.ai_family = AF_UNSPEC;      // don't care IPv4 or IPv6
    hints.ai_socktype = SOCK_STREAM;  // TCP stream sockets
    hints.ai_flags = AI_PASSIVE;      // fill in my IP for me

    getaddrinfo(NULL, port_string.c_str(), &hints, &res);
    server_sock_ = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
    bind(server_sock_, res->ai_addr, res->ai_addrlen);
    listen(server_sock_, 1);
    conn_sock_ = accept(server_sock_, (struct sockaddr *)&their_addr, &addr_size);

    sendln("0");
    sendln(std::to_string(L));
    sendln(std::to_string(N));
    sendln(std::to_string(N_p));
    for (int i = 0; i < N_p_; ++i) {
        sendln(std::to_string(mass[i]));
    }

    freeaddrinfo(res); // free the linked-list
}

Server::~Server() {
}

void Server::marshal(const double *position) {
}

void Server::sendln(const std::string& s) {
    const std::string line = s + "\n";
    send(conn_sock_, line.c_str(), line.length(), 0);
}
