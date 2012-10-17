#ifndef NETWORK_H_
#define NETWORK_H_

#include <winsock2.h>
#include <ws2tcpip.h>
#include <iostream>
#include <vector>
#include <string>

struct fd_set;

void w_init();
int w_close(int sockfd);
void w_cleanup();
std::vector<std::string> w_getnameinfo_vec(int sockfd, int flags=0);
int w_getaddrinfo(const char *node, const char *service,
			  const struct addrinfo *hints, struct addrinfo **res);
int w_socket(int domain, int type, int protocol);
int w_connect(int sockfd, struct sockaddr *serv_addr, int addr_len);
int w_bind(int sockfd, struct sockaddr *my_addr, int addr_len);
int w_listen(int sockfd, int backlog);
int w_accept(int sockfd, struct sockaddr *addr, int *addr_len);
int w_send(int sockfd, char *buf, size_t len, int flags);
int w_sendall(int sockfd, char *buf, unsigned long *len);
int w_recv(int sockfd, char *buf, size_t len, int flags);
int w_recvall(int sockfd, char *buf, unsigned long *len);
int w_select(int numfds, fd_set *readfds, fd_set *writefds,
		   fd_set *exceptfds, struct timeval *timeout);
int w_memcpy_s(void *dest, size_t number_of_elements, const void *src, size_t count);
//addrinfo* w_bind_first_avl(addrinfo *servinfo)
void w_print_servinfo(struct addrinfo *res, std::ostream &fout);

#endif /* NETWORK_H_ */